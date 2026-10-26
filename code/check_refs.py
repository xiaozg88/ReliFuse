# -*- coding: utf-8 -*-
import os
import re

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAIN = os.path.join(BASE, "paper", "v51_ReliFuse_paper.md")
SUPP = os.path.join(BASE, "paper", "v51_ReliFuse_Supplementary.md")


def roman_to_int(s):
    vals = {"I": 1, "V": 5, "X": 10}
    total = 0
    prev = 0
    for ch in reversed(s):
        v = vals[ch]
        if v < prev:
            total -= v
        else:
            total += v
            prev = v
    return total


def int_to_roman(n):
    # enough for 1..20
    table = [(10, "X"), (9, "IX"), (5, "V"), (4, "IV"), (1, "I")]
    out = ""
    for num, sym in table:
        while n >= num:
            out += sym
            n -= num
    return out


def check_main():
    s = open(MAIN, encoding="utf-8").read()
    lines = s.splitlines()
    caps = []
    cap_line_nums = []
    for i, line in enumerate(lines, 1):
        m = re.match(r"^Table ([IVX]+)\.\s", line.strip())
        if m:
            caps.append(m.group(1))
            cap_line_nums.append(i)

    refs = {}
    for i, line in enumerate(lines, 1):
        if re.match(r"^Table [IVX]+\.\s", line.strip()):
            continue
        for m in re.finditer(r"Table ([IVX]+)(?![A-Za-z])", line):
            refs.setdefault(m.group(1), []).append(i)

    print("MAIN table caption order:", caps)
    print("MAIN table first caption line:", cap_line_nums)
    print("MAIN table referenced numbers:", sorted(set(refs), key=roman_to_int))
    for num in caps:
        if num not in refs:
            print("  !! table caption not referenced:", num)
    for num in refs:
        if num not in caps:
            print("  !! reference to nonexistent table:", num)

    # figures
    fig_caps = []
    for line in lines:
        m = re.search(r"Fig\.\s*(\d+)\.", line)
        if m and line.strip().startswith("!"):
            fig_caps.append(int(m.group(1)))
    fig_refs = set()
    for line in lines:
        if line.strip().startswith("!["):
            continue
        for m in re.finditer(r"Fig\.\s*(\d+)([ab])?", line):
            fig_refs.add(int(m.group(1)))
    print("MAIN figure caption numbers:", fig_caps)
    print("MAIN figure referenced numbers:", sorted(fig_refs))
    for n in fig_caps:
        if n not in fig_refs:
            print("  !! figure caption not referenced:", n)
    for n in fig_refs:
        if n not in fig_caps:
            print("  !! reference to nonexistent figure:", n)

    sup_refs = sorted(set(int(x) for x in re.findall(r"Supplementary Table (\d+)", s)))
    print("MAIN references to Supplementary Table:", sup_refs)


def check_supp():
    s = open(SUPP, encoding="utf-8").read()
    lines = s.splitlines()
    t_head = []
    f_head = []
    for line in lines:
        m = re.match(r"^## \u8865\u5145\u8868 (\d+)\.", line.strip())
        if m:
            t_head.append(int(m.group(1)))
        m = re.match(r"^### \u8865\u5145\u56fe (\d+)\.", line.strip())
        if m:
            f_head.append(int(m.group(1)))
    print("SUPP table headings:", t_head)
    print("SUPP figure headings:", f_head)
    print("SUPP table missing:", sorted(set(range(0, max(t_head) + 1)) - set(t_head)))
    print("SUPP figure missing:", sorted(set(range(1, max(f_head) + 1)) - set(f_head)))

    t_refs = set(int(x) for x in re.findall(r"\u8865\u5145\u8868 (\d+)", s))
    f_refs = set(int(x) for x in re.findall(r"\u8865\u5145\u56fe (\d+)", s))
    # headings themselves also match; subtract headings
    print("SUPP missing table targets referenced:", sorted(t_refs - set(t_head)))
    print("SUPP missing figure targets referenced:", sorted(f_refs - set(f_head)))

    main_tabs = set(int_to_roman(i) for i in range(1, 16))
    # refs to main tables in supp: 正文表 X
    zt = set(re.findall(r"\u6b63\u6587\u8868 ([IVX]+)", s))
    print("SUPP references to main tables:", sorted(zt, key=roman_to_int))
    print("SUPP main-table refs not in I-XV:", sorted(zt - main_tabs, key=roman_to_int))
    zf = set(int(x) for x in re.findall(r"\u6b63\u6587\u56fe (\d+)", s))
    print("SUPP references to main figures:", sorted(zf))
    print("SUPP main-figure refs outside 1-6:", sorted(zf - set(range(1, 7))))


if __name__ == "__main__":
    check_main()
    check_supp()
