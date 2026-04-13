#!/usr/bin/env python3
"""Generate SVG vector graphics illustrating the principle behind the FFT transform.

Three diagrams are produced:
  1. fft_signal_decomposition.svg  — time-domain: a composite waveform built
     from three pure sine waves, each component shown alongside the sum.
  2. fft_spectrum.svg              — frequency-domain: the DFT magnitude spectrum
     of the composite signal, showing the three peaks that correspond to the
     three frequency components.
  3. fft_butterfly.svg             — the 8-point radix-2 Cooley–Tukey butterfly
     computation graph, illustrating how FFT recursively splits the DFT.

Run from the repository root:
    python3 tools/generate_fft_svgs.py
Output files are written to docs/.
"""

import cmath
import math
import os


# ---------------------------------------------------------------------------
# Tiny SVG builder
# ---------------------------------------------------------------------------

class SVG:
    """Minimal SVG document builder."""

    def __init__(self, width: int, height: int, title: str = ""):
        self._w = width
        self._h = height
        self._title = title
        self._els: list[str] = []

    # -- primitives ----------------------------------------------------------

    def line(self, x1, y1, x2, y2, stroke="#333", sw=1.5, **extra):
        attrs = f'x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" ' \
                f'stroke="{stroke}" stroke-width="{sw}"'
        for k, v in extra.items():
            attrs += f' {k.replace("_", "-")}="{v}"'
        self._els.append(f"<line {attrs}/>")

    def polyline(self, points: list[tuple[float, float]], stroke="#333",
                 sw=1.5, fill="none", **extra):
        pts = " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
        attrs = f'points="{pts}" stroke="{stroke}" stroke-width="{sw}" fill="{fill}"'
        for k, v in extra.items():
            attrs += f' {k.replace("_", "-")}="{v}"'
        self._els.append(f"<polyline {attrs}/>")

    def rect(self, x, y, w, h, fill="#4c8", stroke="none", sw=1, rx=0, **extra):
        attrs = (f'x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" '
                 f'fill="{fill}" stroke="{stroke}" stroke-width="{sw}" rx="{rx}"')
        for k, v in extra.items():
            attrs += f' {k.replace("_", "-")}="{v}"'
        self._els.append(f"<rect {attrs}/>")

    def circle(self, cx, cy, r, fill="#fff", stroke="#333", sw=1.5):
        self._els.append(
            f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{r:.2f}" '
            f'fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')

    def text(self, x, y, content, *, font_size=13, anchor="middle",
             fill="#222", bold=False, italic=False):
        style = f"font-size:{font_size}px;text-anchor:{anchor};fill:{fill};"
        if bold:
            style += "font-weight:bold;"
        if italic:
            style += "font-style:italic;"
        self._els.append(
            f'<text x="{x:.2f}" y="{y:.2f}" style="{style}">{content}</text>')

    def arrow(self, x1, y1, x2, y2, stroke="#555", sw=1.5, marker="url(#arr)"):
        self._els.append(
            f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" '
            f'stroke="{stroke}" stroke-width="{sw}" marker-end="{marker}"/>')

    def group(self, content: str, **attrs):
        attr_str = " ".join(f'{k.replace("_", "-")}="{v}"' for k, v in attrs.items())
        self._els.append(f"<g {attr_str}>{content}</g>")

    # -- output --------------------------------------------------------------

    def render(self) -> str:
        defs = (
            '<defs>'
            '<marker id="arr" markerWidth="8" markerHeight="8" '
            'refX="6" refY="3" orient="auto">'
            '<path d="M0,0 L0,6 L8,3 z" fill="#555"/>'
            '</marker>'
            '</defs>'
        )
        body = "\n  ".join(self._els)
        title_el = f"<title>{self._title}</title>" if self._title else ""
        return (
            '<?xml version="1.0" encoding="UTF-8"?>\n'
            f'<svg xmlns="http://www.w3.org/2000/svg" '
            f'width="{self._w}" height="{self._h}" '
            f'viewBox="0 0 {self._w} {self._h}" '
            f'font-family="sans-serif">\n'
            f'  {title_el}\n'
            f'  {defs}\n'
            f'  {body}\n'
            '</svg>\n'
        )

    def save(self, path: str):
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(self.render())
        print(f"  wrote {path}")


# ---------------------------------------------------------------------------
# 1. Signal decomposition diagram
# ---------------------------------------------------------------------------

def make_signal_decomposition(path: str):
    W, H = 820, 560
    svg = SVG(W, H, "FFT: Signal Decomposition")

    # white background
    svg.rect(0, 0, W, H, fill="#ffffff", stroke="none")

    # -- title
    svg.text(W / 2, 32, "Signal Decomposition — building a composite waveform",
             font_size=17, bold=True, fill="#111")

    # three component frequencies (normalised periods relative to plot width)
    components = [
        (1,   0.55, "#e05", "f₁ = 1 Hz"),
        (2.5, 0.35, "#07b", "f₂ = 2.5 Hz"),
        (5,   0.25, "#2a7", "f₃ = 5 Hz"),
    ]

    N_PTS = 400
    PAD_L, PAD_R = 70, 30
    row_h = 110       # height per panel
    gap = 14          # gap between panels
    panel_tops = [60, 60 + row_h + gap, 60 + 2 * (row_h + gap),
                  60 + 3 * (row_h + gap)]

    plot_w = W - PAD_L - PAD_R

    def draw_panel(top, amp_scale, signal_fn, color, label, axis_label):
        mid = top + row_h / 2
        bot = top + row_h
        # bounding box (light fill)
        svg.rect(PAD_L, top, plot_w, row_h, fill="#f7f8fa",
                 stroke="#ccd", sw=1, rx=4)
        # zero line
        svg.line(PAD_L, mid, PAD_L + plot_w, mid,
                 stroke="#bbb", sw=1, stroke_dasharray="4,3")
        # waveform
        pts = []
        for i in range(N_PTS + 1):
            t = i / N_PTS           # 0..1  ≡ one "screen"
            y_val = signal_fn(t)
            x = PAD_L + t * plot_w
            y = mid - y_val * amp_scale
            pts.append((x, y))
        svg.polyline(pts, stroke=color, sw=2.2)
        # axes
        svg.arrow(PAD_L, mid, PAD_L + plot_w + 16, mid)                  # time
        svg.arrow(PAD_L, bot, PAD_L, top - 6)                            # amp
        svg.text(PAD_L + plot_w + 22, mid + 4, "t", font_size=12,
                 anchor="start", fill="#555", italic=True)
        svg.text(PAD_L - 8, top - 10, axis_label, font_size=11,
                 anchor="end", fill="#555")
        # label
        svg.text(PAD_L + 10, top + 16, label, font_size=12,
                 anchor="start", fill=color, bold=True)

    # individual component panels
    for idx, (freq, amp, color, label) in enumerate(components):
        draw_panel(
            panel_tops[idx],
            amp_scale=row_h * 0.38,
            signal_fn=lambda t, f=freq, a=amp: a * math.sin(2 * math.pi * f * t),
            color=color,
            label=label,
            axis_label="A",
        )

    # composite sum panel
    def composite(t):
        return sum(a * math.sin(2 * math.pi * f * t)
                   for f, a, *_ in components)

    draw_panel(
        panel_tops[3],
        amp_scale=row_h * 0.32,
        signal_fn=composite,
        color="#333",
        label="x(t) = A₁·sin(2πf₁t) + A₂·sin(2πf₂t) + A₃·sin(2πf₃t)",
        axis_label="A",
    )

    # "+" and "=" connectors on the left margin
    for i in range(len(components) - 1):
        cy = panel_tops[i] + row_h + gap / 2
        svg.text(PAD_L / 2, cy + 5, "+", font_size=18, fill="#888", bold=True)
    cy = panel_tops[2] + row_h + gap / 2
    svg.text(PAD_L / 2, cy + 5, "=", font_size=18, fill="#888", bold=True)

    svg.save(path)


# ---------------------------------------------------------------------------
# 2. Frequency-domain spectrum diagram
# ---------------------------------------------------------------------------

def _dft(x: list[float]) -> list[complex]:
    N = len(x)
    return [sum(x[n] * cmath.exp(-2j * math.pi * k * n / N)
                for n in range(N)) / N
            for k in range(N)]


def make_spectrum(path: str):
    W, H = 760, 420
    svg = SVG(W, H, "FFT: Frequency-Domain Spectrum")

    svg.rect(0, 0, W, H, fill="#ffffff", stroke="none")
    svg.text(W / 2, 32, "Frequency-Domain Magnitude Spectrum  |X(f)|",
             font_size=17, bold=True, fill="#111")

    N = 128
    freqs_amps = [(1, 0.55), (2.5, 0.35), (5, 0.25)]
    signal = [sum(a * math.sin(2 * math.pi * f * n / N)
                  for f, a in freqs_amps)
              for n in range(N)]

    spectrum = _dft(signal)
    mags = [abs(c) for c in spectrum[: N // 2]]
    max_mag = max(mags) if max(mags) > 0 else 1

    PAD_L, PAD_R, PAD_T, PAD_B = 70, 30, 60, 60
    plot_w = W - PAD_L - PAD_R
    plot_h = H - PAD_T - PAD_B
    ox = PAD_L
    oy = PAD_T + plot_h

    # grid & axes
    svg.rect(ox, PAD_T, plot_w, plot_h, fill="#f7f8fa", stroke="#dde", sw=1, rx=4)
    for g in [0.25, 0.5, 0.75, 1.0]:
        gy = oy - g * plot_h
        svg.line(ox, gy, ox + plot_w, gy, stroke="#dde", sw=1)
        svg.text(ox - 6, gy + 4, f"{g:.2f}", font_size=10, anchor="end", fill="#888")

    svg.arrow(ox, oy, ox + plot_w + 16, oy)
    svg.arrow(ox, oy, ox, PAD_T - 10)
    svg.text(ox + plot_w + 24, oy + 4, "f (Hz)", font_size=12,
             anchor="start", fill="#555", italic=True)
    svg.text(ox - 8, PAD_T - 14, "|X(f)|", font_size=12,
             anchor="end", fill="#555")

    # bars
    bar_w = plot_w / (N // 2) * 0.7
    colors = {1: "#e05", 3: "#07b", 6: "#2a7"}   # approximate DFT bins
    freq_labels = {1: "f₁", 3: "f₂", 6: "f₃"}

    # find dominant bins
    dominant = sorted(range(len(mags)), key=lambda i: -mags[i])[:3]

    for k, mag in enumerate(mags):
        bar_h = (mag / max_mag) * plot_h * 0.88
        bx = ox + k * plot_w / (N // 2)
        color = "#4c8e9f"
        if k in dominant:
            rank = dominant.index(k)
            color = ["#e05", "#07b", "#2a7"][rank]
        if bar_h > 1:
            svg.rect(bx - bar_w / 2, oy - bar_h, bar_w, bar_h,
                     fill=color, stroke="none", rx=1)

    # x-axis tick labels (show 0..N//4)
    for k in range(0, N // 4 + 1, N // 16):
        tx = ox + k * plot_w / (N // 2)
        svg.line(tx, oy, tx, oy + 4, stroke="#999", sw=1)
        svg.text(tx, oy + 15, str(k), font_size=10, fill="#666")

    # annotate dominant bins
    for rank, k in enumerate(dominant):
        bx = ox + k * plot_w / (N // 2)
        bar_h = (mags[k] / max_mag) * plot_h * 0.88
        label = [" f₁", " f₂", " f₃"][rank]
        color = ["#e05", "#07b", "#2a7"][rank]
        svg.text(bx, oy - bar_h - 8, label, font_size=12,
                 anchor="middle", fill=color, bold=True)

    # caption
    svg.text(W / 2, H - 12,
             "The three spectral peaks reveal the three hidden frequency components in x(t).",
             font_size=11, fill="#666")

    svg.save(path)


# ---------------------------------------------------------------------------
# 3. Butterfly (radix-2 Cooley–Tukey) diagram
# ---------------------------------------------------------------------------

def make_butterfly(path: str):
    """Draw an 8-point radix-2 DIT FFT butterfly graph."""
    W, H = 800, 520
    svg = SVG(W, H, "FFT: 8-point Radix-2 Butterfly")

    svg.rect(0, 0, W, H, fill="#ffffff", stroke="none")
    svg.text(W / 2, 28, "8-point Radix-2 Cooley–Tukey FFT — Butterfly Graph",
             font_size=17, bold=True, fill="#111")

    N = 8
    # bit-reversal permutation of input indices
    def bit_rev(x, bits):
        r = 0
        for _ in range(bits):
            r = (r << 1) | (x & 1)
            x >>= 1
        return r

    bits = 3   # log2(8)
    order = [bit_rev(i, bits) for i in range(N)]  # [0,4,2,6,1,5,3,7]

    # layout
    n_stages = bits + 1   # input + 3 butterfly stages = 4 columns total
    PAD_L, PAD_R = 80, 80
    PAD_T, PAD_B = 60, 80
    col_w = (W - PAD_L - PAD_R) / (n_stages - 1)
    row_h = (H - PAD_T - PAD_B) / (N - 1)

    def node_xy(stage, row):
        return PAD_L + stage * col_w, PAD_T + row * row_h

    # ---- draw twiddle-factor labels per stage ----
    stage_labels = [
        "Input\n(bit-reversed)",
        "Stage 1\nW⁰",
        "Stage 2\nW⁰, W²",
        "Stage 3\nW⁰…W³",
        "Output\nX[k]",
    ]

    # ---- draw butterfly connections ----
    # Radix-2 DIT: at stage s (0-indexed), group size = 2^(s+1),
    # butterfly span = 2^s.
    twiddle_colors = ["#b06020", "#1a6eb5", "#1a8040"]

    for s in range(bits):
        group = 1 << (s + 1)    # elements per group
        half = group >> 1
        for g in range(N // group):
            for k in range(half):
                top_row = g * group + k
                bot_row = top_row + half
                x1, y1 = node_xy(s, top_row)
                x2, y2 = node_xy(s, bot_row)
                xo1, yo1 = node_xy(s + 1, top_row)
                xo2, yo2 = node_xy(s + 1, bot_row)
                col = twiddle_colors[s % len(twiddle_colors)]
                # top input → top output (straight)
                svg.line(x1, y1, xo1, yo1, stroke=col, sw=1.4)
                # top input → bottom output (cross)
                svg.line(x1, y1, xo2, yo2, stroke=col, sw=1.4,
                         stroke_dasharray="5,3")
                # bottom input → top output (cross, negated)
                svg.line(x2, y2, xo1, yo1, stroke=col, sw=1.4,
                         stroke_dasharray="5,3")
                # bottom input → bottom output (straight, negated)
                svg.line(x2, y2, xo2, yo2, stroke=col, sw=1.4)

                # twiddle W^k label on the cross arm (top → bottom)
                mx = (x1 + xo2) / 2 + 4
                my = (y1 + yo2) / 2
                svg.text(mx, my - 2, f"W^{k}", font_size=9,
                         anchor="start", fill=col)

    # ---- draw nodes ----
    for s in range(n_stages):
        for r in range(N):
            x, y = node_xy(s, r)
            svg.circle(x, y, 7, fill="#fff", stroke="#334", sw=1.5)

    # ---- input labels ----
    for r in range(N):
        x, y = node_xy(0, r)
        svg.text(x - 14, y + 4, f"x[{order[r]}]", font_size=11,
                 anchor="end", fill="#333")

    # ---- output labels ----
    for r in range(N):
        x, y = node_xy(n_stages - 1, r)
        svg.text(x + 14, y + 4, f"X[{r}]", font_size=11,
                 anchor="start", fill="#333")

    # ---- stage column headers ----
    header_y = PAD_T - 28
    # 4 columns: Input, Stage 1, Stage 2, Stage 3 (= output after all 3 butterfly stages)
    stage_short = ["Input", "Stage 1", "Stage 2", "Stage 3"]
    for s, label in enumerate(stage_short):
        x, _ = node_xy(s, 0)
        svg.text(x, header_y, label, font_size=11,
                 fill="#555", bold=(s in (0, n_stages - 1)))

    # ---- butterfly legend ----
    lx, ly = PAD_L + 10, H - 55
    svg.text(lx, ly, "butterfly unit:", font_size=11, anchor="start", fill="#444")
    svg.line(lx + 110, ly - 4, lx + 160, ly - 4, stroke="#555", sw=1.4)
    svg.text(lx + 165, ly, "a + W^k·b", font_size=10, anchor="start", fill="#555")
    svg.line(lx + 110, ly + 10, lx + 160, ly + 10, stroke="#555", sw=1.4,
             stroke_dasharray="5,3")
    svg.text(lx + 165, ly + 14, "a − W^k·b", font_size=10, anchor="start", fill="#555")

    svg.text(W / 2, H - 12,
             "W^k = e^(−2πik/N)  (twiddle factor).  "
             "Dashed lines carry a sign flip; solid lines do not.",
             font_size=11, fill="#666")

    svg.save(path)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    out_dir = os.path.join(os.path.dirname(__file__), "..", "docs")
    print("Generating FFT diagrams …")
    make_signal_decomposition(os.path.join(out_dir, "fft_signal_decomposition.svg"))
    make_spectrum(os.path.join(out_dir, "fft_spectrum.svg"))
    make_butterfly(os.path.join(out_dir, "fft_butterfly.svg"))
    print("Done.")
