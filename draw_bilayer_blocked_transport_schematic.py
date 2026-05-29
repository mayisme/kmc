#!/usr/bin/env python3
from __future__ import annotations

import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

OUT = Path('/home/xiao/文档/KMC/comparison_al2o3_zno_vs_10nm_al2o3_envelope/bilayer_6nm_zno_blocked_transport_schematic_4k.png')

W, H = 3000, 4200
M = 210


def font(size: int, bold: bool = False):
    base = '/usr/share/fonts/truetype/dejavu/DejaVuSans'
    return ImageFont.truetype(f'{base}-Bold.ttf' if bold else f'{base}.ttf', size=size)


F_TITLE = font(86, True)
F_HEAD = font(58, True)
F_LABEL = font(44, False)
F_BOLD = font(46, True)
F_SMALL = font(36, False)
F_TINY = font(30, False)


def text_center(draw, xy, text, fnt, fill='black', spacing=8):
    x, y = xy
    lines = text.split('\n')
    bbs = [draw.textbbox((0, 0), line, font=fnt) for line in lines]
    total_h = sum(bb[3] - bb[1] for bb in bbs) + spacing * (len(lines) - 1)
    cy = y - total_h / 2
    for line, bb in zip(lines, bbs):
        w = bb[2] - bb[0]
        h = bb[3] - bb[1]
        draw.text((x - w / 2, cy), line, font=fnt, fill=fill)
        cy += h + spacing


def arrow(draw, start, end, fill='#1267c7', width=9, head=30):
    draw.line([start, end], fill=fill, width=width)
    x1, y1 = start
    x2, y2 = end
    a = math.atan2(y2 - y1, x2 - x1)
    pts = [
        (x2, y2),
        (x2 - head * math.cos(a - 0.48), y2 - head * math.sin(a - 0.48)),
        (x2 - head * math.cos(a + 0.48), y2 - head * math.sin(a + 0.48)),
    ]
    draw.polygon(pts, fill=fill)


def water(draw, x, y, r=18):
    draw.ellipse((x - r, y - r, x + r, y + r), fill='#0c7fe8', outline='#004d9c', width=4)
    draw.ellipse((x - r * 0.35, y - r * 0.35, x + r * 0.15, y + r * 0.15), fill='#92d9ff')


def blob(draw, cx, cy, rx, ry, fill, outline='#53607d', width=4, seed=0):
    pts = []
    for i in range(64):
        a = 2 * math.pi * i / 64
        noise = 1 + 0.13 * math.sin(3 * a + seed) + 0.08 * math.cos(7 * a + seed * 0.7)
        pts.append((cx + rx * noise * math.cos(a), cy + ry * noise * math.sin(a)))
    draw.polygon(pts, fill=fill, outline=outline)
    if width > 1:
        draw.line(pts + [pts[0]], fill=outline, width=width)


def wavy_polygon(x0, y0, x1, y1, amp_top=16, amp_bot=14, n=120):
    top = []
    bot = []
    for i in range(n):
        t = i / (n - 1)
        x = x0 + (x1 - x0) * t
        top.append((x, y0 + amp_top * math.sin(2 * math.pi * (6 * t + 0.15)) + 5 * math.sin(2 * math.pi * 17 * t)))
        bot.append((x, y1 + amp_bot * math.sin(2 * math.pi * (5 * t + 0.4))))
    return top + list(reversed(bot))


def draw_scale_bar(draw, x, y, width, label, color):
    draw.line((x, y, x + width, y), fill=color, width=8)
    draw.line((x, y - 20, x, y + 20), fill=color, width=8)
    draw.line((x + width, y - 20, x + width, y + 20), fill=color, width=8)
    draw.text((x + width / 2, y + 34), label, font=F_SMALL, fill=color, anchor='ma')


def cross_block(draw, x, y, r=46):
    draw.ellipse((x - r, y - r, x + r, y + r), fill='white', outline='#cf0000', width=8)
    draw.line((x - r * 0.55, y - r * 0.55, x + r * 0.55, y + r * 0.55), fill='#cf0000', width=9)
    draw.line((x - r * 0.55, y + r * 0.55, x + r * 0.55, y - r * 0.55), fill='#cf0000', width=9)


def draw_path(draw, pts, color='#e3342f', width=9, dot=True):
    draw.line(pts, fill=color, width=width, joint='curve')
    for a, b in zip(pts[:-1], pts[1:]):
        arrow(draw, a, b, fill=color, width=0, head=0)
    if dot:
        for i, (x, y) in enumerate(pts):
            draw.ellipse((x - 14, y - 14, x + 14, y + 14), fill='white', outline=color, width=5)


def main():
    img = Image.new('RGB', (W, H), 'white')
    d = ImageDraw.Draw(img)

    d.text((W / 2, 88), 'Moisture transport in 4.5 nm Al2O3 / 6 nm ZnO bilayer', font=F_TITLE, fill='#101820', anchor='ma')
    d.text((W / 2, 178), '6 nm ZnO is thinner than a typical ~10 nm grain: draw one truncated lateral grain layer, not stacked grains',
           font=F_LABEL, fill='#34495e', anchor='ma')

    body_x0, body_x1 = 620, 2440
    al_y0, al_y1 = 610, 1730
    zno_y0, zno_y1 = al_y1, 2530

    d.text((W / 2, 360), 'Water vapor (H2O)', font=F_HEAD, fill='#0057c2', anchor='ma')
    for i in range(34):
        x = body_x0 + 80 + (i * 137) % (body_x1 - body_x0 - 160)
        y = 430 + 55 * (i % 4)
        water(d, x, y, 17)

    # Al2O3 layer
    al_poly = wavy_polygon(body_x0, al_y0, body_x1, al_y1)
    d.polygon(al_poly, fill='#c7cbe8', outline='#6974a8')
    d.line(al_poly + [al_poly[0]], fill='#6b73a6', width=5)

    # porous/free-volume features in Al2O3
    pores = [
        (760, 760, 95, 60), (1020, 885, 120, 78), (1360, 720, 82, 55),
        (1660, 940, 135, 88), (2100, 790, 105, 70), (2280, 1210, 120, 83),
        (840, 1300, 110, 75), (1215, 1460, 78, 54), (1545, 1325, 95, 66),
        (1900, 1490, 125, 85), (1370, 1035, 60, 45),
    ]
    for idx, p in enumerate(pores):
        blob(d, *p, fill='#8790c5', outline='#eef1ff', width=5, seed=idx + 1)
    # blocked pinhole-like defect that terminates before interface
    blob(d, 1475, 620, 48, 32, fill='#b3f3f1', outline='#315b6d', width=5, seed=20)
    blob(d, 1508, 850, 58, 45, fill='#b3f3f1', outline='#315b6d', width=5, seed=21)
    blob(d, 1450, 1115, 50, 42, fill='#b3f3f1', outline='#315b6d', width=5, seed=22)
    cross_block(d, 1500, 1255, 42)

    # interface
    d.line((body_x0, al_y1, body_x1, al_y1), fill='#2e7d67', width=7)

    # ZnO single truncated grain layer
    d.rectangle((body_x0, zno_y0, body_x1, zno_y1), fill='#c9efe2', outline='#2c6f5e', width=5)
    boundaries = [
        [(820, zno_y0), (880, 1840), (840, 2020), (910, zno_y1)],
        [(1160, zno_y0), (1100, 1880), (1215, 2110), (1160, zno_y1)],
        [(1510, zno_y0), (1590, 1910), (1500, 2210), (1540, zno_y1)],
        [(1880, zno_y0), (1815, 1880), (1935, 2050), (1865, zno_y1)],
        [(2220, zno_y0), (2265, 1900), (2170, 2180), (2230, zno_y1)],
    ]
    for pts in boundaries:
        d.line(pts, fill='#2c7a6b', width=13, joint='curve')
        d.line(pts, fill='#0f4d42', width=4, joint='curve')
    # Subtle large-grain interiors. These are clipped by the 6 nm film surfaces,
    # so they should read as one laterally spread grain row, not stacked grains.
    for i, x in enumerate([720, 1010, 1350, 1720, 2070, 2360]):
        blob(d, x, 2120, 145, 220, fill='#d9f6ed', outline='#a8d9ca', width=3, seed=60 + i)
    # Redraw grain boundaries over the pale interiors.
    for pts in boundaries:
        d.line(pts, fill='#2c7a6b', width=13, joint='curve')
        d.line(pts, fill='#0f4d42', width=4, joint='curve')

    # blocked grain boundary pathway
    path1 = [(1475, 390), (1510, 620), (1470, 810), (1520, 1000), (1480, 1160), (1505, 1250)]
    d.line(path1, fill='#e3342f', width=9, joint='curve')
    for x, y in path1:
        d.ellipse((x - 15, y - 15, x + 15, y + 15), fill='white', outline='#e3342f', width=5)
    d.text((1585, 1210), 'blocked inside\nAl2O3 defect\n(no continuous pore)', font=F_SMALL, fill='#c1121f')
    arrow(d, (1580, 1205), (1538, 1250), fill='#c1121f', width=5, head=24)

    path2 = [(1040, 390), (1110, 620), (1125, 845), (1090, 1110), (1165, 1430), (1260, 1710),
             (1320, 1840), (1260, 1990), (1345, 2140)]
    d.line(path2, fill='#e3342f', width=9, joint='curve')
    for x, y in path2:
        d.ellipse((x - 15, y - 15, x + 15, y + 15), fill='white', outline='#e3342f', width=5)
    cross_block(d, 1345, 2140, 46)
    d.text((1415, 2070), 'ZnO grain boundary\nis sealed / dead-end;\npath does not reach\nreceiver side', font=F_SMALL, fill='#c1121f')
    arrow(d, (1415, 2075), (1365, 2140), fill='#c1121f', width=5, head=24)

    # a third water molecule stopped at interface mismatch
    path3 = [(2010, 395), (2070, 620), (2040, 840), (2090, 1090), (2035, 1370), (2110, 1640)]
    d.line(path3, fill='#e3342f', width=8, joint='curve')
    for x, y in path3:
        d.ellipse((x - 13, y - 13, x + 13, y + 13), fill='white', outline='#e3342f', width=4)
    cross_block(d, 2110, 1660, 38)
    d.text((2150, 1550), 'blocked at\nAl2O3/ZnO\ninterface mismatch', font=F_TINY, fill='#c1121f')

    # labels left
    d.text((170, 925), 'Amorphous\nAl2O3', font=F_HEAD, fill='#2647a8')
    d.line((300, 1100, 540, 1100), fill='#2647a8', width=5)
    d.text((110, 1175), '4.5 nm layer;\ndisordered free volume\nand local defects,\nbut no guaranteed\nthrough-pore', font=F_SMALL, fill='#2c3e50')

    d.text((190, 2090), 'Polycrystalline\nZnO', font=F_HEAD, fill='#087a59')
    d.line((300, 2295, 540, 2295), fill='#087a59', width=5)
    d.text((110, 2365), '6 nm layer;\nthinner than typical\n~10 nm lateral grains;\ntransport only along\ngrain boundaries', font=F_SMALL, fill='#2c3e50')

    # labels right
    d.text((2490, 905), 'Pores / free volume\nin amorphous Al2O3', font=F_SMALL, fill='#30395c')
    d.line((2365, 895, 2470, 895), fill='#30395c', width=4)
    d.text((2490, 1370), 'Terminated defect\nblocks H2O path', font=F_SMALL, fill='#30395c')
    d.line((2240, 1360, 2470, 1360), fill='#30395c', width=4)
    d.text((2490, 1945), 'One-row truncated\nZnO grains', font=F_SMALL, fill='#0f4d42')
    d.line((2210, 1940, 2470, 1940), fill='#0f4d42', width=4)
    d.text((2490, 2245), 'Grain boundary\nblocked / dead-end', font=F_SMALL, fill='#0f4d42')
    d.line((1940, 2235, 2470, 2235), fill='#0f4d42', width=4)

    # scale bars
    draw_scale_bar(d, 2550, al_y0, 190, '4.5 nm Al2O3', '#2647a8')
    draw_scale_bar(d, 2550, zno_y0, 250, '6 nm ZnO', '#087a59')
    d.text((2475, 2580), 'Schematic not to scale laterally:\nZnO grain width ~10 nm,\nfilm thickness only 6 nm.', font=F_TINY, fill='#2c3e50')

    # Receiver side, no arrival
    d.text((W / 2, 2755), 'Receiver side: no transmitted H2O in blocked-path scenario', font=F_HEAD, fill='#0057c2', anchor='ma')
    d.line((body_x0, 2675, body_x1, 2675), fill='#d7e6ff', width=4)

    # legend and interpretation
    legend_y = 2970
    d.rounded_rectangle((M, legend_y, W - M, H - 190), radius=36, fill='#fbfbfb', outline='#b7b7b7', width=3)
    water(d, M + 70, legend_y + 85, 18)
    d.text((M + 115, legend_y + 58), '= H2O molecule', font=F_SMALL, fill='black')
    d.line((M + 460, legend_y + 85, M + 560, legend_y + 85), fill='#e3342f', width=9)
    d.ellipse((M + 503, legend_y + 72, M + 529, legend_y + 98), fill='white', outline='#e3342f', width=4)
    d.text((M + 590, legend_y + 58), '= attempted KMC path', font=F_SMALL, fill='black')
    cross_block(d, M + 1030, legend_y + 85, 28)
    d.text((M + 1075, legend_y + 58), '= blocked / no transmission', font=F_SMALL, fill='black')
    blob(d, M + 70, legend_y + 180, 30, 22, fill='#b3f3f1', outline='#315b6d', width=4, seed=40)
    d.text((M + 115, legend_y + 154), '= defect / pinhole / free-volume pocket', font=F_SMALL, fill='black')
    d.line((M + 920, legend_y + 180, M + 1045, legend_y + 180), fill='#0f4d42', width=10)
    d.text((M + 1080, legend_y + 154), '= ZnO grain boundary', font=F_SMALL, fill='black')

    interp_lines = [
        'Key idea: because the ZnO film is only 6 nm thick, a ~10 nm ZnO grain is cut by the film surfaces.',
        'The barrier benefit is not many vertical grain layers, but the lack of aligned through-defects:',
        'Al2O3 defects can terminate, and ZnO grain-boundary pathways can be blocked or dead-ended before reaching the receiver side.',
    ]
    yy = legend_y + 285
    for line in interp_lines:
        d.text((M + 70, yy), line, font=F_SMALL, fill='#243447')
        yy += 54

    OUT.parent.mkdir(parents=True, exist_ok=True)
    img.save(OUT)
    print(OUT)


if __name__ == '__main__':
    main()
