#!/usr/bin/env python3
"""Mechanical half of the `main` → `programs::<name>::run(argv, out)` conversion (3.10).

Moves a binary's source into `src/programs/<category>/` — the same category the binary sits in
under `src/bin/` — inserts the `Program` record, turns `main` into `run`, and redirects every
`println!`/`print!` into the `out` sink.

Two shapes of binary, told apart by the file itself:

  * `fn main()` delegating to a `fn run() -> Result<(), String>` that already parses with
    `Arguments` — converted whole: `USAGE` and the accepted keys are read out of the source, so
    neither has to be passed in.
  * the compute in `main()` — converted as far as it goes; the hand-rolled argument loop becomes
    `Arguments` by hand afterwards, because only the author knows which keys are multi-valued,
    which are flags and what the defaults are.

Usage: convert_program.py <crate> <bin-relative-path> <name> <status> [known,keys] [<<<'usage']
"""

import pathlib
import re
import subprocess
import sys


def convert_prints(src: str) -> str:
    """`println!(a, b)` → `writeln!(out, a, b).map_err(io_err)?`, matching parens.

    `eprintln!`/`eprint!` are left alone: stderr is not part of the `out` contract.
    """
    for macro, repl in (("println!", "writeln!"), ("print!", "write!")):
        i = 0
        while True:
            m = re.search(r"(?<![a-z_])" + re.escape(macro) + r"\(", src[i:])
            if not m:
                break
            start = i + m.start()
            open_paren = i + m.end() - 1
            depth, j, in_str, esc = 1, open_paren + 1, False, False
            while j < len(src) and depth:
                c = src[j]
                if esc:
                    esc = False
                elif c == "\\":
                    esc = True
                elif c == '"':
                    in_str = not in_str
                elif not in_str:
                    depth += (c == "(") - (c == ")")
                j += 1
            args = src[open_paren + 1 : j - 1]
            sink = "out" if not args.strip() else "out, "
            src = src[:start] + f"{repl}({sink}{args}).map_err(io_err)?" + src[j:]
            i = start + 1
    return src


# `box` is a directory under `src/bin/` and a keyword: its module is `r#box`.
KEYWORDS = {"box", "move", "type", "match", "loop", "ref", "const", "fn", "mod", "use", "where"}


def ident(category: str) -> str:
    """The module path for a category directory."""
    return f"r#{category}" if category in KEYWORDS else category


def register_in_category(root: pathlib.Path, category: str, name: str) -> list[str]:
    """Add `pub mod <name>;` to the category's mod.rs. Returns every program it now holds."""
    mod = root / "src/programs" / category / "mod.rs"
    text = mod.read_text() if mod.exists() else f"//! The {category} programs.\n\n"
    names = sorted(set(re.findall(r"^pub mod (\w+);", text, flags=re.M)) | {name})
    body = "".join(f"pub mod {x};\n" for x in names)
    if re.search(r"^pub mod \w+;\n", text, flags=re.M):
        text = re.sub(r"(^pub mod \w+;\n)+", body, text, count=1, flags=re.M)
    else:
        text = text.rstrip("\n") + "\n\n" + body
    mod.write_text(text)
    return names


def insert_program_entry(m: str, name: str, siblings: list[str], category: str) -> str:
    """Put `&<name>::PROGRAM,` among its category's entries, keeping that section sorted."""
    lines = m.split("\n")
    entry = f"    &{name}::PROGRAM,"
    at = [i for i, l in enumerate(lines) if l.strip() in {f"&{s}::PROGRAM," for s in siblings}]
    if not at:  # a category with no entry yet: a new section before the closing bracket
        close = next(i for i, l in enumerate(lines) if l.startswith("];"))
        lines[close:close] = ([f"    // {category}"] if category else []) + [entry]
        return "\n".join(lines)
    after = [i for i in at if lines[i].strip() > entry.strip()]
    lines.insert(after[0] if after else at[-1] + 1, entry)
    return "\n".join(lines)


def add_reexport(m: str, path: str, name: str) -> str:
    """Keep the program reachable as `programs::<name>`, not only through its category."""
    group = re.search(rf"pub use {re.escape(path)}::\{{([^}}]*)\}};", m)
    if group:
        names = sorted(set(re.findall(r"\w+", group.group(1))) | {name})
        return m[: group.start()] + f"pub use {path}::{{{', '.join(names)}}};" + m[group.end() :]
    single = re.search(rf"pub use {re.escape(path)}::(\w+);", m)
    if single:
        names = sorted({single.group(1), name})
        return m[: single.start()] + f"pub use {path}::{{{', '.join(names)}}};" + m[single.end() :]
    uses = list(re.finditer(r"^pub use .*;\n", m, flags=re.M))
    return splice(m, f"pub use {path}::{name};\n", uses[-1].end() if uses else None)


def splice(m: str, line: str, at: int | None = None) -> str:
    """Insert `line` at `at`, or after the `pub mod` block — after the imports if there is none."""
    if at is None:
        block = re.search(r"(^pub mod .*;\n)+", m, flags=re.M)
        head = re.search(r"^use .*;\n", m, flags=re.M)
        at = block.end() if block else head.end() + 1
    return m[:at] + line + m[at:]


def drop_import(text: str, item: str) -> str:
    """Remove `item` from the file's `use` statements; it went out with the code that used it."""
    text = re.sub(rf"^use [\w:]+::{item};\n", "", text, flags=re.M)
    text = re.sub(rf"^\s*{item},\n", "", text, flags=re.M)  # its own line in a braced list
    text = re.sub(rf"\{{{item}, ", "{", text)
    text = re.sub(rf", {item}\}}", "}", text)
    return re.sub(r"\{(\w+)\};", r"\1;", text)  # a lone survivor needs no braces


def add_imports(text: str, category: str) -> str:
    """The three things a `run` needs that a `main` did not, each only if it is missing."""
    wanted = []
    if "io_err" not in text:
        wanted.append(f"use {'crate::programs' if category else 'super'}::io_err;")
    if not re.search(r"\bStatus\b", text):
        wanted.append("use newen_io::args::{Program, Status};")
    if not re.search(r"^use std::io::(\{[^}]*\bWrite\b|Write;)", text, flags=re.M):
        wanted.append("use std::io::Write;")
    if not wanted:
        return text
    lines = text.split("\n")
    last_use = max((i for i, l in enumerate(lines) if l.startswith("use ")), default=-1)
    lines[last_use + 1 : last_use + 1] = wanted
    return "\n".join(lines)


def convert_delegating(text: str, name: str, status: str, category: str) -> str:
    """A `main` that only forwards to `run`: drop it and let `run` take the argv itself."""
    text = re.sub(r"^fn main\(\) \{\n(?:.*\n)*?^\}\n", "", text, count=1, flags=re.M)
    text = text.replace(
        "fn run() -> Result<(), String> {",
        "pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {",
        1,
    )

    # the keys are already in the source, as the array `Arguments::parse` is called with
    call = re.search(r"Arguments::parse\(\s*&\[(.*?)\],\s*USAGE,?\s*\)", text, flags=re.S)
    if not call:
        raise SystemExit(f"{name}: no `Arguments::parse(&[...], USAGE)` to take the keys from")
    keys = re.findall(r'"([^"]+)"', call.group(1))
    text = (
        text[: call.start()]
        + "Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)"
        + text[call.end() :]
    )

    text = re.sub(r"^const USAGE:", "pub const USAGE:", text, flags=re.M)
    text = drop_import(text, "fail") if "fail(" not in text else text
    text = add_imports(text, category)

    record = (
        f"\npub const KNOWN: &[&str] = &[{', '.join(f'\"{k}\"' for k in keys)}];\n\n"
        "pub const PROGRAM: Program = Program {\n"
        f'    name: "{name}",\n    known: KNOWN,\n    usage: USAGE,\n'
        f"    status: Status::{status},\n    run,\n}};"
    )
    # after the USAGE literal, which runs to the line ending it
    lines = text.split("\n")
    start = next(i for i, l in enumerate(lines) if l.startswith("pub const USAGE:"))
    end = next(i for i in range(start, len(lines)) if lines[i].rstrip().endswith('";'))
    lines.insert(end + 1, record)
    return "\n".join(lines)


def main() -> int:
    crate, bin_path, name, status = sys.argv[1:5]
    known = sys.argv[5] if len(sys.argv) > 5 else ""  # read from the source when the shape allows
    usage = "" if sys.stdin.isatty() else sys.stdin.read().rstrip("\n")

    root = pathlib.Path(crate)
    lib = root.name.replace("-", "_")
    src_path = root / "src/bin" / bin_path
    parent = pathlib.PurePosixPath(bin_path).parent.as_posix()
    category = "" if parent == "." else parent  # `src/bin/` is flat in some crates
    dst_dir = root / "src/programs" / category
    dst_dir.mkdir(parents=True, exist_ok=True)
    dst_path = dst_dir / f"{name}.rs"
    subprocess.run(["git", "mv", str(src_path), str(dst_path)], check=True)

    text = dst_path.read_text()

    # the body lives in the library now, where the crate cannot name itself
    text = text.replace(f"{lib}::", "crate::")

    # drop the process-owning imports; the program never touches any of them
    text = re.sub(r"^use std::env;\n", "", text, flags=re.M)
    text = re.sub(r"^use std::process;\n", "", text, flags=re.M)
    text = drop_import(text, "gromos_args")

    if "fn run() -> Result<(), String> {" in text:
        text = convert_delegating(text, name, status, category)
        text = convert_prints(text)
        dst_path.write_text(text)
        write_shim(src_path, lib, name)
        register(root, category, name)
        print(f"{name}: converted, {len(text.splitlines())} lines")
        return 0

    keys = ", ".join(f'"{k}"' for k in known.split(",") if k)
    header = (
        f"use {'crate::programs' if category else 'super'}::io_err;\n"
        "use newen_io::args::{Arguments, Program, Status};\n"
        "use std::io::Write;\n"
    )
    literal = '"' + usage.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n") + '"'
    record = (
        f"\npub const USAGE: &str = {literal};\n\n"
        + f"pub const KNOWN: &[&str] = &[{keys}];\n\n"
        + "pub const PROGRAM: Program = Program {\n"
        + f'    name: "{name}",\n    known: KNOWN,\n    usage: USAGE,\n'
        + f"    status: Status::{status},\n    run,\n}};\n"
    )

    # insert the header after the last top-level `use`, the record after it; a `use` can span
    # lines, so the insertion point is the line that *ends* the last one
    lines = text.split("\n")
    depth, last_use = 0, -1
    for i, l in enumerate(lines):
        if l.startswith("use ") or depth:
            depth += l.count("{") - l.count("}")
            if not depth:
                last_use = i
    if last_use < 0:
        doc_end = max((i for i, l in enumerate(lines) if l.startswith("//!")), default=-1)
        lines.insert(doc_end + 1, "\n" + header.rstrip() + "\n" + record)
    else:
        lines.insert(last_use + 1, header.rstrip() + "\n" + record)
    text = "\n".join(lines)

    text = text.replace(
        "fn main() {",
        "pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {",
        1,
    )
    text = convert_prints(text)
    text = close_with_ok(text)
    dst_path.write_text(text)

    write_shim(src_path, lib, name)
    register(root, category, name)

    print(f"{name}: moved, {len(text.splitlines())} lines — now convert its argument parsing")
    return 0


def close_with_ok(text: str) -> str:
    """`main` returned nothing; `run` returns `Result`, so the body needs a tail `Ok(())`."""
    text = text.rstrip("\n") + "\n"
    end = text.rindex("}\n")
    return text if text[:end].rstrip().endswith("Ok(())") else text[:end] + "    Ok(())\n}\n"


def write_shim(src_path: pathlib.Path, lib: str, name: str) -> None:
    """What is left of the binary: the one place allowed to own the process."""
    src_path.write_text(
        f"//! `{name}` — the body is `{lib}::programs::{name}`.\n\n"
        f"fn main() {{\n"
        f"    newen_io::run_main(&{lib}::programs::{name}::PROGRAM)\n"
        f"}}\n"
    )


def register(root: pathlib.Path, category: str, name: str) -> None:
    """A `pub mod` in its category, the flat re-export and a PROGRAMS entry at the root."""
    mod = root / "src/programs/mod.rs"
    m = mod.read_text()
    if f"&{name}::PROGRAM," not in m:
        if category:
            siblings = register_in_category(root, category, name)
            path = ident(category)
            if not re.search(rf"^pub mod {re.escape(path)};$", m, flags=re.M):
                m = splice(m, f"pub mod {path};\n")
            m = add_reexport(m, path, name)
        else:
            siblings = sorted(re.findall(r"^pub mod (\w+);", m, flags=re.M) + [name])
            m = re.sub(
                r"(^pub mod \w+;\n)+",
                "".join(f"pub mod {x};\n" for x in siblings),
                m,
                count=1,
                flags=re.M,
            )
        m = insert_program_entry(m, name, siblings, category)
        mod.write_text(m)


if __name__ == "__main__":
    sys.exit(main())
