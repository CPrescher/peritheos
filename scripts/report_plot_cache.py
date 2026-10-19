"""Integrity-checked plot inputs, with bounded cross-platform replay roundoff.

Snapshots retain full numerical precision. The tolerance is only for deciding
whether a saved figure is current; model evaluation and source data are unchanged.
"""

from __future__ import annotations

import base64
import binascii
import hashlib
import json
import math
import re
import zlib
from xml.etree import ElementTree

import numpy as np


def same_inputs(left, right):
    if isinstance(left, dict) and isinstance(right, dict):
        return left.keys() == right.keys() and all(
            same_inputs(left[k], right[k]) for k in left
        )
    if isinstance(left, (list, tuple)) and isinstance(right, (list, tuple)):
        if len(left) != len(right):
            return False
        # Only convert flat numeric vectors. Older supported NumPy versions
        # warn instead of rejecting ragged lists of differently sized panels.
        if all(
            isinstance(v, (int, float)) and not isinstance(v, bool)
            for values in (left, right)
            for v in values
        ):
            return bool(np.allclose(left, right, rtol=1e-9, atol=1e-9, equal_nan=False))
        return all(same_inputs(a, b) for a, b in zip(left, right))
    if isinstance(left, bool) or isinstance(right, bool):
        return type(left) is type(right) and left == right
    if isinstance(left, (float, int)) and isinstance(right, (float, int)):
        return math.isclose(left, right, rel_tol=1e-9, abs_tol=1e-9)
    return type(left) is type(right) and left == right


def inputs_json(inputs):
    return json.dumps(
        inputs, default=lambda a: a.tolist(), allow_nan=False, sort_keys=True
    )


def finish_svg(payload, inputs, code_signature):
    payload = b"\n".join(line.rstrip() for line in payload.splitlines()) + b"\n"
    snapshot = base64.b64encode(zlib.compress(inputs_json(inputs).encode(), level=9))
    payload = payload.replace(
        b"?>\n", b"?>\n<!-- peritheos-plot-inputs " + snapshot + b" -->\n", 1
    )
    marker = f"<!-- peritheos-curves code={code_signature} body={hashlib.sha256(payload).hexdigest()} -->\n".encode()
    return payload.replace(b"?>\n", b"?>\n" + marker, 1)


def current_svg(payload, inputs, code_signature):
    marker = re.search(
        rb"<!-- peritheos-curves code=([a-f0-9]{64}) body=([a-f0-9]{64}) -->\n",
        payload,
    )
    if not marker or marker[1].decode() != code_signature:
        return False
    body = payload[: marker.start()] + payload[marker.end() :]
    if hashlib.sha256(body).hexdigest() != marker[2].decode():
        return False
    snapshot = re.search(rb"<!-- peritheos-plot-inputs ([A-Za-z0-9+/=]+) -->", body)
    if not snapshot:
        return False
    try:
        saved = json.loads(
            zlib.decompress(base64.b64decode(snapshot[1], validate=True))
        )
        return ElementTree.fromstring(
            payload
        ).tag == "{http://www.w3.org/2000/svg}svg" and same_inputs(
            saved, json.loads(inputs_json(inputs))
        )
    except (ValueError, zlib.error, binascii.Error, ElementTree.ParseError):
        return False
