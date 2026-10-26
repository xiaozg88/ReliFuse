# -*- coding: utf-8 -*-
"""Regenerate the JBHI LaTeX from the v51 markdown source."""
import io
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
MD = os.path.join(ROOT, "paper", "v51_ReliFuse_paper.md")
TEX_DIR = os.path.join(ROOT, "paper", "JBHI_template", "IEEE-TJ-color-latex-template")
TEX = os.path.join(TEX_DIR, "ReliFuse.tex")

# ---------------------------------------------------------------------------
# equation LaTeX forms (keys match the markdown [[EQ:...:n]] numbers)
EQ_LATEX = {
    1: r"\kappa_{1}(i,j)=\frac{1+\cos(x_{i},x_{j})}{2}",
    2: r"\kappa_{2}(i,j)=\frac{1}{1+|i-j|}",
    3: r"\kappa_{3}(i,j)=\frac{\min(d_{i},d_{j})}{\max(d_{i},d_{j})}",
    4: r"\kappa_{4}(i,j)=\exp\!\left(-\frac{\lVert x_{i}-x_{j}\rVert^{2}}{d}\right)",
    5: r"g_{e}^{(l)}=\beta\sum_{k=1}^{4}\alpha_{l,b,k}\,\kappa_{k}(e)+(1-\beta)",
    6: r"z_{mc}=\frac{z_{m}}{T_{m}}, \quad m\in\{1,2\}",
    7: r"p_{m}=\sigma(z_{mc})=\frac{1}{1+\exp(-z_{mc})}",
    8: r"H(p)=-[p\ln p+(1-p)\ln(1-p)]",
    9: r"\begin{aligned}\mathrm{agreement}&=p_{1}p_{2}+(1-p_{1})(1-p_{2}),\\ \mathrm{disagreement}&=|p_{1}-p_{2}|\end{aligned}",
    10: r"x_{\mathrm{gate}}=[z,p,H,\mathrm{agreement},\mathrm{disagreement},s_{\mathrm{struct}},\rho]",
    11: r"g=\sigma(\mathrm{MLP}_{g}(h)),\quad b=\mathrm{MLP}_{b}(h)",
    12: r"z_{\mathrm{fuse}}=z_{1c}+g(z_{2c}-z_{1c})+b",
}

# figure label -> (pdf filename, width style)
FIG_MAP = {
    "fig1_architecture": ("fig1_architecture.pdf", "0.85\\textwidth"),
    "fig2_protocol": ("fig2_protocol.pdf", "0.85\\textwidth"),
    "fig_two_axis_complementarity": ("fig_two_axis_complementarity.png", "linewidth"),
    "fig3_main_roc_pr": ("fig3_main_roc_pr.pdf", "0.85\\textwidth"),
    "fig7_fusion_anatomy": ("fig7_fusion_anatomy.pdf", "0.9\\textwidth"),
    "fig6_mechanism_kernels": ("fig6_mechanism_kernels.pdf", "0.9\\textwidth"),
}

SECTION_LABEL = {
    "Introduction": "sec:introduction",
    "Related Work": "sec:related",
    "Materials and Methods": "sec:methods",
    "Experimental Results": "sec:results",
    "Discussion": "sec:discussion",
    "Conclusion": "sec:conclusion",
}


def convert_citations(text):
    # [1], [2], [3] -> \cite{b1,b2,b3};  [1] -> \cite{b1}
    def repl(m):
        nums = re.findall(r"\[(\d+)\]", m.group(0))
        # Preserve a single trailing space if the original citation group was
        # followed by whitespace, so the next word does not run into the cite.
        tail = " " if m.group(0).endswith((" ", "\t", "\n")) else ""
        return "\\cite{" + ",".join("b" + n for n in nums) + "}" + tail
    return re.sub(r"(?:\[\d+\]\s*)+", repl, text)


def to_latex_text(text):
    text = convert_citations(text)
    # Escape the LaTeX comment character before any other processing so that
    # strings such as "95% CI" or "10–30%" do not truncate the rest of a line.
    text = text.replace("%", "\\%")
    # literal underscore subscripts (d_i, z_s, s_struct, f_low, p_m, ...) -> inline math
    text = re.sub(r"\b([A-Za-z])_([A-Za-z0-9]+)\b", r"$\1_{\2}$", text)
    # multi-letter / digit subscripted identifiers
    text = text.replace("margin_m", "$\\mathrm{margin}_{m}$")
    text = text.replace("margin_1", "$\\mathrm{margin}_{1}$").replace("margin_2", "$\\mathrm{margin}_{2}$")
    text = text.replace("z1c", "$z_{1c}$").replace("z2c", "$z_{2c}$")
    text = text.replace("p1", "$p_{1}$").replace("p2", "$p_{2}$")
    # unicode Greek / math symbols -> inline math
    text = text.replace("α", "$\\alpha$").replace("β", "$\\beta$")
    text = text.replace("Δ", "$\\Delta$").replace("σ", "$\\sigma$")
    text = text.replace("κ", "$\\kappa$").replace("ρ", "$\\rho$")
    text = text.replace("→", "$\\to$").replace("≈", "$\\approx$").replace("≥", "$\\ge$")
    text = text.replace("×", "$\\times$").replace("±", "$\\pm$")
    # less/greater-than and absolute-value bars (OT1 renders them as i/¿/em-dash)
    text = text.replace("<", "$<$").replace(">", "$>$").replace("|", "$|$")
    # superscripts
    text = text.replace("²", "${}^2$").replace("¹", "${}^1$")
    # dashes
    text = text.replace("−", "$-$").replace("–", "--").replace("—", "---")
    # Angstrom
    text = text.replace("Å", "{\\AA}")
    # accented letters (rare, in author names)
    text = text.replace("ö", "\\\"{o}").replace("č", "\\v{c}").replace("ć", "\\'{c}").replace("ò", "\\`{o}")
    # Markdown inline formatting used by the v51 sources.
    text = re.sub(r"`([^`]+)`", r"\\texttt{\1}", text)
    text = re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", text)
    text = re.sub(r"\*([^*]+?)\*", r"\\emph{\1}", text)
    return text


def markdown_table_to_latex(rows):
    # rows: list of list of cells; skip separator row
    rows = [r for r in rows if not all(re.fullmatch(r":?-{2,}:?", c.strip()) for c in r)]
    if not rows:
        return ""
    header = rows[0]
    ncol = len(header)
    body = rows[1:]
    colspec = "l" + "c" * (ncol - 1)
    out = ["\\begin{tabular}{" + colspec + "}", "\\hline"]
    out.append(" & ".join(to_latex_text(c.strip()) for c in header) + " \\\\")
    out.append("\\hline")
    for r in body:
        cells = [to_latex_text(c.strip()) for c in r]
        # pad to ncol
        while len(cells) < ncol:
            cells.append("")
        out.append(" & ".join(cells) + " \\\\")
    out.append("\\hline")
    out.append("\\end{tabular}")
    return "\n".join(out)


def build_body():
    md = io.open(MD, encoding="utf-8").read()
    lines = md.splitlines()
    out = []
    # skip everything before "## I. Introduction"
    start = next(i for i, l in enumerate(lines) if l.strip().startswith("## I. Introduction"))
    # stop before "## References"
    stop = next(i for i, l in enumerate(lines) if l.strip().startswith("## References"))

    i = start
    pending_caption = None
    while i < stop:
        line = lines[i]
        s = line.strip()
        if not s:
            i += 1
            continue
        # section
        m = re.match(r"^## (.+)$", s)
        if m:
            title = m.group(1)
            # strip leading Roman numeral "I. " etc.
            title = re.sub(r"^[IVX]+\.\s*", "", title)
            if title in SECTION_LABEL:
                out.append("\\section{" + title + "}\n\\label{" + SECTION_LABEL[title] + "}")
            else:
                out.append("\\section*{" + title + "}")
            i += 1
            continue
        # subsection
        m = re.match(r"^### ([A-Z])\.\s*(.+)$", s)
        if m:
            out.append("\\subsection{" + m.group(2) + "}")
            i += 1
            continue
        # equation
        m = re.match(r"^\[\[EQ:.*?:(\d+)\]\]$", s)
        if m:
            k = int(m.group(1))
            out.append("\\begin{equation}\n" + EQ_LATEX[k] + "\n\\label{eq:" + str(k) + "}\n\\end{equation}")
            i += 1
            continue
        # figure
        m = re.match(r"^!\[(.*)\]\((.*)\)$", s)
        if m:
            raw_caption = m.group(1)
            # figure number comes from the original caption before cleaning
            fn = re.search(r"Fig\.\s*(\d+)", raw_caption)
            num = fn.group(1) if fn else "?"
            # IEEEtran already prints "Fig. N."; avoid duplicating the number.
            caption = re.sub(r"^Fig\.\s*\d+\.\s*", "", raw_caption)
            path = m.group(2)
            base = os.path.basename(path).replace(".png", "").replace(".pdf", "")
            fname, width = FIG_MAP.get(base, (os.path.basename(path), "linewidth"))
            if width == "linewidth":
                block = ("\\begin{figure}[!t]\n\\centering\n"
                         "\\includegraphics[width=\\linewidth]{../../../figures/" + fname + "}\n"
                         "\\caption{" + to_latex_text(caption) + "}\n"
                         "\\label{fig:" + num + "}\n\\end{figure}")
            else:
                block = ("\\begin{figure*}[!t]\n"
                         "\\centerline{\\includegraphics[width=" + width + "]{../../../figures/" + fname + "}}\n"
                         "\\caption{" + to_latex_text(caption) + "}\n"
                         "\\label{fig:" + num + "}\n\\end{figure*}")
            out.append(block)
            i += 1
            continue
        # table caption line "Table I. ..."
        if re.match(r"^Table [IVX]+\.\s", s) or re.match(r"^Table \d+\.\s", s):
            pending_caption = s
            i += 1
            continue
        # table (markdown pipe)
        if s.startswith("|"):
            tbl = []
            while i < stop and lines[i].strip().startswith("|"):
                cells = [c.strip() for c in lines[i].strip().strip("|").split("|")]
                tbl.append(cells)
                i += 1
            tabular = markdown_table_to_latex(tbl)
            ncol = max(len(r) for r in tbl if not all(re.fullmatch(r":?-{2,}:?", c.strip()) for c in r))
            tab_block = tabular
            if ncol >= 99:
                tab_block = "\\resizebox{\\textwidth}{!}{" + tabular + "}"
            table_env = "table*" if ncol >= 5 else "table"
            block = "\\begin{" + table_env + "}[!t]\n"
            if pending_caption:
                block += "\\caption{" + to_latex_text(pending_caption) + "}\n"
                tn = re.search(r"Table ([IVX]+|\d+)", pending_caption)
                if tn:
                    block += "\\label{tab:" + tn.group(1).lower() + "}\n"
                pending_caption = None
            block += ("\\centering\n\\setlength{\\tabcolsep}{3pt}\\small\n" + tab_block + "\n\\end{" + table_env + "}")
            out.append(block)
            continue
        # bullet
        if s.startswith("- ") or s.startswith("* "):
            items = []
            while i < stop and lines[i].strip().startswith(("- ", "* ")):
                items.append("\\item " + to_latex_text(lines[i].strip()[2:]))
                i += 1
            out.append("\\begin{itemize}\n" + "\n".join(items) + "\n\\end{itemize}")
            continue
        # normal paragraph
        out.append(to_latex_text(s))
        i += 1

    return "\n\n".join(out)


def main():
    cur = io.open(TEX, encoding="utf-8").read()
    header = cur[:cur.find("\\end{IEEEkeywords}") + len("\\end{IEEEkeywords}")]
    refs_start = cur.find("\\begin{thebibliography}")
    refs = cur[refs_start:] if refs_start != -1 else "\\begin{thebibliography}{33}\\end{thebibliography}\n\\end{document}"
    body = build_body()
    full = header + "\n\n" + body + "\n\n" + refs
    io.open(TEX, "w", encoding="utf-8").write(full)
    print("wrote", TEX, "chars:", len(full))


if __name__ == "__main__":
    main()
