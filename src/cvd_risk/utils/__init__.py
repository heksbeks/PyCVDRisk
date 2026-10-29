"""Utility functions for CVD risk calculator."""

import math
from pathlib import Path

def round_half_up(x, ndigits=None):
	'''implements standard rounding'''
	if not math.isfinite(x):
		return x
	d = -0.5 if x < 0 else 0.5
	if ndigits is None:
		return int(math.floor(x + d))
	m = 10**ndigits
	return math.floor(x*m + d)/m

def getLatestFolder(root, pattern):
	dirs = [f for f in Path(root).glob(pattern) if f.is_dir()]
	if not dirs:
		return None
	return str(max(dirs, key = lambda f: f.stat().st_birthtime))
