"""
extract_coastline.py
====================
Extracts the Western Australian coastline used in Figure 7 from the geo-maps
``countries-land-100m`` dataset and clips it to the study bounding box.

Data provenance
---------------
Coastline polygons are from the geo-maps ``countries-land-100m`` package
(Primarosa, 2020; https://github.com/simonepri/geo-maps), which derives its
geometry from Natural Earth 1:10 m physical vectors and OpenStreetMap. In the
manuscript workflow the polygons are handled with the sf package in R
(Pebesma, 2018); this script provides an equivalent dependency-free extraction
so that the figure can be reproduced from Python alone.

Usage
-----
    npm install @geo-maps/countries-land-100m
    python extract_coastline.py \
        node_modules/@geo-maps/countries-land-100m/map.geo.json \
        ../results/coast_wa.json

The 45 KB ``coast_wa.json`` produced here is shipped with the repository so the
figures render without the 121 MB world file.
"""
import json
import sys

# study bounding box (a little wider than the visible map so only true coast shows)
XMIN, XMAX, YMIN, YMAX = 114.0, 116.6, -34.6, -29.4


def clip_rect(ring):
    """Sutherland-Hodgman clip of a ring ([[x, y], ...]) to the bounding box."""
    def clip_edge(pts, inside, intersect):
        out, n = [], len(pts)
        if n == 0:
            return out
        for i in range(n):
            a, b = pts[i], pts[(i + 1) % n]
            ina, inb = inside(a), inside(b)
            if ina and inb:
                out.append(b)
            elif ina and not inb:
                out.append(intersect(a, b))
            elif not ina and inb:
                out.append(intersect(a, b))
                out.append(b)
        return out

    def ix(a, b, xv):
        t = (xv - a[0]) / (b[0] - a[0]) if b[0] != a[0] else 0.0
        return [xv, a[1] + t * (b[1] - a[1])]

    def iy(a, b, yv):
        t = (yv - a[1]) / (b[1] - a[1]) if b[1] != a[1] else 0.0
        return [a[0] + t * (b[0] - a[0]), yv]

    p = ring
    p = clip_edge(p, lambda q: q[0] >= XMIN, lambda a, b: ix(a, b, XMIN))
    p = clip_edge(p, lambda q: q[0] <= XMAX, lambda a, b: ix(a, b, XMAX))
    p = clip_edge(p, lambda q: q[1] >= YMIN, lambda a, b: iy(a, b, YMIN))
    p = clip_edge(p, lambda q: q[1] <= YMAX, lambda a, b: iy(a, b, YMAX))
    return p


def main(src, dst, a3="AUS"):
    gj = json.load(open(src))
    feats = gj["features"] if gj.get("type") == "FeatureCollection" else [gj]
    target = [f for f in feats if (f.get("properties") or {}).get("A3") == a3]
    g = target[0]["geometry"]
    polys = g["coordinates"] if g["type"] == "MultiPolygon" else [g["coordinates"]]

    clipped = []
    for poly in polys:
        for ring in poly:
            xs = [c[0] for c in ring]; ys = [c[1] for c in ring]
            if max(xs) < XMIN or min(xs) > XMAX or max(ys) < YMIN or min(ys) > YMAX:
                continue
            cr = clip_rect([[c[0], c[1]] for c in ring])
            if len(cr) >= 3:
                clipped.append([[round(x, 5), round(y, 5)] for x, y in cr])
    clipped.sort(key=len, reverse=True)
    json.dump(clipped, open(dst, "w"))
    print(f"wrote {dst}: {len(clipped)} land polygons (largest {len(clipped[0])} pts)")


if __name__ == "__main__":
    src = sys.argv[1] if len(sys.argv) > 1 else \
        "node_modules/@geo-maps/countries-land-100m/map.geo.json"
    dst = sys.argv[2] if len(sys.argv) > 2 else "../results/coast_wa.json"
    main(src, dst)
