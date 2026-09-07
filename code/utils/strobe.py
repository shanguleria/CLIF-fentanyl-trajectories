"""CONSORT/STROBE cohort flow: text and figure.

FLOW rows are dicts with step, n_before, n_excluded, reason, n_after. Exclusions
are derived from consecutive counts rather than asserted, so the diagram cannot
disagree with the numbers it was built from.
"""
from __future__ import annotations

from pathlib import Path

INK = "#1a1a1a"
BOX = "#f4f6f8"
EDGE = "#5b6b7a"
EXCL = "#fdf3f2"
EXCL_EDGE = "#b4736c"


def render_text(flow: list[dict]) -> str:
    w = max(len(r["step"]) for r in flow) + 2
    out = ["CONSORT / STROBE cohort flow", "=" * (w + 26), ""]
    for i, r in enumerate(flow):
        out.append(f"{r['step']:<{w}} {r['n_after']:>10,}")
        if i + 1 < len(flow) and flow[i + 1]["n_excluded"]:
            nxt = flow[i + 1]
            out += [f"{'':<{w}}     |",
                    f"{'':<{w}}     |-- excluded {nxt['n_excluded']:>8,}  {nxt['reason']}",
                    f"{'':<{w}}     v"]
    return "\n".join(out)


def render_png(flow: list[dict], path: Path, title: str = "") -> Path:
    """Draw the flow. Stages down the left, exclusions branching right.

    Boxes grow to fit wrapped text, so a long step name cannot collide with its
    count and a long reason cannot overflow its box.
    """
    import textwrap

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    STEP_W, EXCL_W = 6.6, 6.2
    LABEL_CHARS, REASON_CHARS = 40, 48
    LINE_H, PAD = 0.30, 0.26
    BX, GAP = 0.35, 0.4
    WIDTH = BX + STEP_W + GAP + EXCL_W + BX

    steps = list(flow)
    lab = [textwrap.wrap(s["step"], LABEL_CHARS) or [""] for s in steps]
    box_h = [len(l) * LINE_H + PAD for l in lab]
    rsn = [textwrap.wrap(s["reason"] or "", REASON_CHARS) or [""] for s in steps]
    exc_h = [(len(r) + 1) * LINE_H + PAD for r in rsn]

    total = box_h[0]
    for i in range(1, len(steps)):
        total += (max(exc_h[i] + 0.5, 0.95) if steps[i]["n_excluded"] else 0.55)
        total += box_h[i]
    height = total + (0.75 if title else 0.35) + 0.25

    fig, ax = plt.subplots(figsize=(WIDTH, max(4.0, height)))
    ax.set_xlim(0, WIDTH)
    ax.set_ylim(0, height)
    ax.axis("off")
    bx = BX

    if title:
        ax.text(bx, height - 0.32, title, fontsize=13, fontweight="bold", color=INK)
    y = height - (0.75 if title else 0.35)

    for i, st in enumerate(steps):
        h = box_h[i]
        y -= h
        ax.add_patch(FancyBboxPatch(
            (bx, y), STEP_W, h, boxstyle="round,pad=0.02,rounding_size=0.06",
            linewidth=1.1, edgecolor=EDGE, facecolor=BOX))
        for k, line in enumerate(lab[i]):
            ax.text(bx + 0.16, y + h - PAD / 2 - LINE_H * (k + 0.5), line,
                    va="center", ha="left", fontsize=10, color=INK)
        ax.text(bx + STEP_W - 0.16, y + h / 2, f"{st['n_after']:,}",
                va="center", ha="right", fontsize=11.5, color=INK, fontweight="bold")

        if i + 1 >= len(steps):
            break
        nxt = steps[i + 1]
        drop = max(exc_h[i + 1] + 0.5, 0.95) if nxt["n_excluded"] else 0.55
        ax.add_patch(FancyArrowPatch(
            (bx + STEP_W / 2, y), (bx + STEP_W / 2, y - drop),
            arrowstyle="-|>", mutation_scale=13, linewidth=1.1, color=EDGE))
        if nxt["n_excluded"]:
            eh = exc_h[i + 1]
            my = y - drop / 2
            ax.add_patch(FancyArrowPatch(
                (bx + STEP_W / 2, my), (bx + STEP_W + GAP - 0.08, my),
                arrowstyle="-|>", mutation_scale=11, linewidth=1.0, color=EXCL_EDGE))
            ex, ey = bx + STEP_W + GAP, my - eh / 2
            ax.add_patch(FancyBboxPatch(
                (ex, ey), EXCL_W, eh,
                boxstyle="round,pad=0.02,rounding_size=0.05",
                linewidth=1.0, edgecolor=EXCL_EDGE, facecolor=EXCL))
            ax.text(ex + 0.16, ey + eh - PAD / 2 - LINE_H * 0.5,
                    f"excluded  {nxt['n_excluded']:,}", va="center", ha="left",
                    fontsize=9.5, color=EXCL_EDGE, fontweight="bold")
            for k, line in enumerate(rsn[i + 1]):
                ax.text(ex + 0.16, ey + eh - PAD / 2 - LINE_H * (k + 1.5), line,
                        va="center", ha="left", fontsize=8.8, color=INK)
        y -= drop

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path
