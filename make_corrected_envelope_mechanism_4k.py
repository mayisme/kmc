#!/usr/bin/env python3
from __future__ import annotations

import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


W, H = 3840, 2160
OUT = Path('/home/xiao/文档/KMC/comparison_al2o3_zno_vs_10nm_al2o3_envelope/al2o3_zno_reference_envelope_mechanism_4k_corrected.png')


def font(size: int, bold: bool = False):
    base = '/usr/share/fonts/truetype/dejavu/DejaVuSans'
    path = f'{base}-Bold.ttf' if bold else f'{base}.ttf'
    return ImageFont.truetype(path, size=size)


F_TITLE = font(58, True)
F_PANEL = font(76, True)
F_HEAD = font(44, True)
F_TEXT = font(36, False)
F_SMALL = font(28, False)
F_BOLD = font(36, True)
F_MATH = font(34, True)


def multiline_center(draw, xy, text, fill, fnt, spacing=8):
    x, y = xy
    lines = text.split('\n')
    bbs = [draw.textbbox((0, 0), line, font=fnt) for line in lines]
    widths = [bb[2] - bb[0] for bb in bbs]
    heights = [bb[3] - bb[1] for bb in bbs]
    total_h = sum(heights) + spacing * (len(lines) - 1)
    cy = y - total_h / 2
    for line, width, height in zip(lines, widths, heights):
        draw.text((x - width / 2, cy), line, font=fnt, fill=fill)
        cy += height + spacing


def arrow(draw, start, end, fill, width=8):
    draw.line([start, end], fill=fill, width=width)
    x1, y1 = start
    x2, y2 = end
    ang = math.atan2(y2 - y1, x2 - x1)
    size = 24
    pts = [
        (x2, y2),
        (x2 - size * math.cos(ang - 0.45), y2 - size * math.sin(ang - 0.45)),
        (x2 - size * math.cos(ang + 0.45), y2 - size * math.sin(ang + 0.45)),
    ]
    draw.polygon(pts, fill=fill)


def water(draw, x, y, r=15):
    draw.ellipse((x - r, y - r, x + r, y + r), fill='#0b7de3', outline='#064a9b', width=3)
    draw.ellipse((x - r / 3, y - r / 3, x + r / 5, y + r / 5), fill='#7ec8ff')


def defect_blob(draw, cx, cy, rx, ry, fill='#a7f0ee', outline='#063a4a'):
    pts = []
    for i in range(40):
        a = 2 * math.pi * i / 40
        noise = 1 + 0.10 * math.sin(3 * a + cx * 0.01) + 0.06 * math.cos(5 * a)
        pts.append((cx + rx * noise * math.cos(a), cy + ry * noise * math.sin(a)))
    draw.polygon(pts, fill=fill, outline=outline)


def layer_rect(draw, box, fill, outline='#101010', width=4, waviness=False):
    x0, y0, x1, y1 = box
    if not waviness:
        draw.rectangle(box, fill=fill, outline=outline, width=width)
        return
    top = []
    bottom = []
    for i in range(80):
        t = i / 79
        x = x0 + t * (x1 - x0)
        top.append((x, y0 + 8 * math.sin(2 * math.pi * (4 * t + 0.1)) + 5 * math.sin(2 * math.pi * 11 * t)))
        bottom.append((x, y1 + 8 * math.sin(2 * math.pi * (5 * t + 0.35))))
    pts = top + list(reversed(bottom))
    draw.polygon(pts, fill=fill, outline=outline)


def draw_zno_grains(draw, box):
    x0, y0, x1, y1 = box
    draw.rectangle(box, fill='#ffe6a3', outline='#101010', width=4)
    cols, rows = 5, 3
    for r in range(rows):
        for c in range(cols):
            cx = x0 + (c + 0.5 + (0.18 if r % 2 else 0)) * (x1 - x0) / cols
            cy = y0 + (r + 0.5) * (y1 - y0) / rows
            radx = (x1 - x0) / cols * 0.55
            rady = (y1 - y0) / rows * 0.58
            pts = []
            for i in range(6):
                a = math.pi / 6 + i * math.pi / 3
                pts.append((cx + radx * math.cos(a), cy + rady * math.sin(a)))
            draw.line(pts + [pts[0]], fill='#e76f00', width=5)


def metric_box(draw, box, lines, color):
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=20, fill='white', outline=color, width=4)
    y = y0 + 26
    for left, right, good in lines:
        draw.text((x0 + 35, y), left, font=F_MATH, fill=color)
        draw.text((x0 + 178, y), right, font=F_BOLD, fill=('#0a7f29' if good else '#c1121f'))
        y += 50


def main() -> None:
    img = Image.new('RGB', (W, H), 'white')
    d = ImageDraw.Draw(img)

    d.text((W // 2, 45), 'Water-vapor barrier mechanism: Al2O3 reference envelope vs Al2O3/ZnO bilayer',
           font=F_TITLE, fill='#111111', anchor='ma')

    margin_x = 70
    top = 160
    panel_h = 1580
    gap = 24
    panel_w = (W - 2 * margin_x - 3 * gap) // 4
    panels = []
    for i in range(4):
        x0 = margin_x + i * (panel_w + gap)
        x1 = x0 + panel_w
        panels.append((x0, top, x1, top + panel_h))
        d.rectangle((x0, top, x1, top + panel_h), outline='#b0b0b0', width=3)
        d.text((x0 + 28, top + 25), chr(ord('A') + i), font=F_PANEL, fill='black')
        d.text(((x0 + x1) // 2, top + 68), 'Humid side', font=F_HEAD, fill='black', anchor='ma')
        d.text(((x0 + x1) // 2, top + panel_h - 325), 'Dry side', font=F_HEAD, fill='black', anchor='ma')
        for k in range(11):
            water(d, x0 + 120 + (k * 67 + 37 * (k % 3)) % (panel_w - 220), top + 190 + 60 * (k % 3), 15)

    al = '#173b69'
    defect = '#a8f0ef'

    # A high-quality dense Al2O3
    x0, y0, x1, y1 = panels[0]
    layer_rect(d, (x0 + 25, y0 + 470, x1 - 25, y0 + 1045), al, waviness=True)
    for cx, cy in [(x0 + 230, y0 + 510), (x0 + 530, y0 + 595), (x0 + 420, y0 + 835), (x0 + 700, y0 + 520)]:
        defect_blob(d, cx, cy, 23, 18, defect)
    multiline_center(d, ((x0 + x1) / 2, y0 + 710), '10 nm dense\namorphous\nAl2O3', 'white', font(52, True))
    arrow(d, (x0 + 190, y0 + 310), (x0 + 230, y0 + 500), '#0b65c8', 6)
    metric_box(d, (x0 + 160, y0 + 1190, x1 - 160, y0 + 1455), [
        ('Ptrans', 'low', True),
        ('WVTR', 'low', True),
        ('defects', 'non-through', True),
    ], '#0b55c0')

    # B through pinhole
    x0, y0, x1, y1 = panels[1]
    layer_rect(d, (x0 + 25, y0 + 470, x1 - 25, y0 + 1045), al, waviness=True)
    d.rounded_rectangle((x0 + 365, y0 + 455, x0 + 485, y0 + 1085), radius=45, fill=defect, outline='#063a4a', width=4)
    for yy in range(y0 + 520, y0 + 1120, 105):
        water(d, x0 + 425, yy, 17)
        arrow(d, (x0 + 425, yy + 28), (x0 + 425, yy + 82), '#0b65c8', 7)
    multiline_center(d, (x0 + 220, y0 + 710), '10 nm\nAl2O3', 'white', font(52, True))
    d.text((x0 + 510, y0 + 455), 'direct\nthrough-pinhole', font=F_TEXT, fill='black')
    arrow(d, (x0 + 500, y0 + 535), (x0 + 455, y0 + 610), 'black', 5)
    metric_box(d, (x0 + 145, y0 + 1155, x1 - 145, y0 + 1515), [
        ('Ptrans', 'high', False),
        ('FPT', 'short', False),
        ('Lpath', 'short', False),
        ('tau', 'low', False),
        ('D_eff', 'high', False),
        ('WVTR', 'high', False),
    ], '#d00000')

    # C staggered non-through
    x0, y0, x1, y1 = panels[2]
    layer_rect(d, (x0 + 25, y0 + 470, x1 - 25, y0 + 1045), al, waviness=True)
    for cx, cy in [(x0 + 330, y0 + 540), (x0 + 500, y0 + 690), (x0 + 690, y0 + 855), (x0 + 260, y0 + 900), (x0 + 610, y0 + 1010)]:
        defect_blob(d, cx, cy, 60, 45, defect)
    multiline_center(d, (x0 + 205, y0 + 710), '10 nm\nAl2O3', 'white', font(52, True))
    d.text((x0 + 570, y0 + 350), 'staggered\nnon-through\ndefects', font=F_TEXT, fill='black')
    arrow(d, (x0 + 570, y0 + 500), (x0 + 500, y0 + 670), 'black', 5)
    metric_box(d, (x0 + 150, y0 + 1200, x1 - 150, y0 + 1475), [
        ('geometry', 'zero-leakage ideal', True),
        ('Ptrans', '~0', True),
        ('WVTR', '~0', True),
    ], '#137d22')

    # D bilayer
    x0, y0, x1, y1 = panels[3]
    al_box = (x0 + 25, y0 + 410, x1 - 25, y0 + 625)
    zno_box = (x0 + 25, y0 + 625, x1 - 25, y0 + 1045)
    layer_rect(d, al_box, al, waviness=True)
    draw_zno_grains(d, zno_box)
    for cx in [x0 + 170, x0 + 350, x0 + 560]:
        d.rounded_rectangle((cx - 18, y0 + 400, cx + 18, y0 + 635), radius=15, fill=defect, outline='#063a4a', width=3)
        arrow(d, (cx, y0 + 305), (cx, y0 + 455), '#0b65c8', 6)
    path = [(x0 + 350, y0 + 635), (x0 + 420, y0 + 715), (x0 + 350, y0 + 800), (x0 + 515, y0 + 900), (x0 + 470, y0 + 1030), (x0 + 500, y0 + 1115)]
    d.line(path, fill='#1455cc', width=9, joint='curve')
    for a, b in zip(path[:-1], path[1:]):
        arrow(d, a, b, '#1455cc', 5)
    multiline_center(d, (x1 - 155, y0 + 500), '4.5 nm\namorphous\nAl2O3', '#0b3f93', F_BOLD)
    multiline_center(d, (x1 - 155, y0 + 815), '6 nm\npolycrystalline\nZnO', '#c45f00', F_BOLD)
    d.text((x0 + 565, y0 + 1000), 'tortuous\ngrain-boundary\nnetwork', font=F_SMALL, fill='black')
    arrow(d, (x0 + 565, y0 + 1005), (x0 + 500, y0 + 930), 'black', 5)
    metric_box(d, (x0 + 145, y0 + 1155, x1 - 145, y0 + 1515), [
        ('FPT', 'up', True),
        ('Lpath', 'up', True),
        ('tau', 'up', True),
        ('D_eff', 'down', True),
        ('WVTR', 'down', True),
    ], '#7b1b84')

    # Legend
    leg_y = 1810
    d.rounded_rectangle((95, leg_y, W - 95, H - 88), radius=26, fill='#fbfbfb', outline='#b5b5b5', width=3)
    water(d, 150, leg_y + 70, 18)
    d.text((190, leg_y + 46), '= water molecule', font=F_TEXT, fill='black')
    defect_blob(d, 150, leg_y + 150, 22, 18, defect)
    d.text((190, leg_y + 126), '= defect / pinhole', font=F_TEXT, fill='black')

    legend = [
        ('Ptrans', '= penetration probability'),
        ('FPT', '= first passage time distribution'),
        ('Lpath', '= water-molecule path length distribution'),
        ('tau', '= tortuosity = Lpath / Lfilm'),
        ('D_eff', '= effective diffusivity'),
        ('WVTR', '= water vapor transmission rate'),
    ]
    lx = 760
    ly = leg_y + 35
    for i, (k, v) in enumerate(legend):
        col = i // 3
        row = i % 3
        x = lx + col * 940
        y = ly + row * 72
        d.text((x, y), k, font=F_MATH, fill='black')
        d.text((x + 145, y), v, font=F_TEXT, fill='black')

    d.text((W // 2, H - 42), 'Performance is an envelope, not a fixed 6185x claim.',
           font=font(42, False), fill='#111111', anchor='ma')

    OUT.parent.mkdir(parents=True, exist_ok=True)
    img.save(OUT)
    print(OUT)


if __name__ == '__main__':
    main()
