#!/usr/bin/env python3
"""Build jsc/main.tex (single file, Springer Nature sn-jnl template) from the
section sources copied verbatim from paper2d/sections.  Mathematical content is
not touched: every section file is inlined byte for byte, except that the
generative-AI statement is moved out of the declarations file into the methods
section (Section 3), as Journal of Scientific Computing's LLM policy requires.

Run from jsc/:  python3 flatten.py   (then build with pdflatex/bibtex, see README.txt)
"""
import pathlib, re

here = pathlib.Path(__file__).resolve().parent
sec = here / "sections"
src_main = (here / "main.src.tex").read_text()
# macro block copied verbatim from the paper2d preamble (kept in jsc/paper2d_main.tex)
p2 = (here / "paper2d_main.tex").read_text()
macros = p2[p2.index("\\newcommand{\\R}"):p2.index("\\begin{document}")].rstrip() + "\n"

# ---- split the declarations file of paper2d into its parts -------------------
decl = (sec / "Z_declarations.tex").read_text()
m = re.search(r"\\subsection\*\{Use of generative AI\}\n(.*?)\n\\subsection\*\{Data and code availability\}", decl, re.S)
ai_text = m.group(1).rstrip()
def part(title, nxt):
    pat = r"\\subsection\*\{" + re.escape(title) + r"\}\n(.*?)" + (r"\n\\subsection\*\{" + re.escape(nxt) + r"\}" if nxt else r"\Z")
    return re.search(pat, decl, re.S).group(1).rstrip()
data_text = part("Data and code availability", "Competing interests")
ci_text = part("Competing interests", "Funding")
fund_text = part("Funding", "Acknowledgements")
ack_text = part("Acknowledgements", None)

ai_section = (
    "\n\\subsection{Use of generative AI}\\label{sec:ai}\n"
    "%% Placed here (the methods section) because the journal's LLM policy asks for LLM use to be\n"
    "%% documented in the Methods section; wording unchanged from paper2d/sections/Z_declarations.tex.\n"
    + ai_text + "\n")

body = []
for name in ["01_intro", "02_related", "03_setting", "04_printed", "05_corrected",
             "06_penalty", "07_strong", "08_further", "09_numerics", "10_conclusions"]:
    t = (sec / f"{name}.tex").read_text().rstrip() + "\n"
    body.append(f"%%%%%%%%%%%%%%%%%%%% {name}.tex (verbatim from paper2d/sections) %%%%%%%%%%%%%%%%%%%%\n" + t)
    if name == "03_setting":
        body.append(ai_section)

ack_block = "\\bmhead{Acknowledgements}\n" + ack_text + "\n"

declarations = r"""\section*{Declarations}

\begin{itemize}
\item \textbf{Funding.} FUNDTEXT
\item \textbf{Competing interests.} CITEXT
\item \textbf{Ethics approval and consent to participate.} Not applicable.
\item \textbf{Consent for publication.} Not applicable.
\item \textbf{Data availability.} DATATEXT
\item \textbf{Materials availability.} Not applicable.
\item \textbf{Code availability.} The finite element code (the Scott--Vogelius--Nitsche solver), the exact rational-arithmetic and interval certificates, and the scripts that produce every table and figure are available at \url{https://github.com/jaideepsaipadhi/zero-gradient-prize} \cite{ZeroGradRepo}; see also \cite{PadhiExt}.
\item \textbf{Author contributions.} J.~S.~Padhi is the sole author and is responsible for the conception of the study, the analysis, the software and computations, and the writing of the manuscript (with the use of generative AI described in Section~\ref{sec:ai}).
\item \textbf{Use of generative AI.} The use of a generative AI model in this work is described in Section~\ref{sec:ai}. The AI model is not an author.
\end{itemize}
"""
declarations = (declarations.replace("FUNDTEXT", fund_text.replace("\n%% VERIFY with the author.", "") +
                                     "\n%% FILL if the author has any funding (grant, fellowship); otherwise keep.")
                .replace("CITEXT", ci_text)
                .replace("DATATEXT", data_text))

out = (src_main.replace("%%MACROS%%", macros).replace("%%BODY%%", "\n".join(body))
               .replace("%%ACK%%", ack_block)
               .replace("%%DECLARATIONS%%", declarations))
(here / "main.tex").write_text(out)
print("wrote main.tex:", len(out.splitlines()), "lines")
