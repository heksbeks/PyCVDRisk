'''
Test Prevent.calculate and Prevent.calculate_batch against preventr
'''

import collections, importlib.util, os, pprint, subprocess, sys, random
from pathlib import Path

script_dir = Path(__file__).parent.resolve()

packages = ['pandas', 'numpy', 'rpy2']
missing = [pkg for pkg in packages if importlib.util.find_spec(pkg.replace('-', '_')) is None]
if missing:
	print(f'Installing missing packages {missing}...')
	subprocess.check_call([sys.executable, '-m', 'pip', 'install'] + missing)
	print('Finished installing packages!\n')

import pandas as pd

sys.path.insert(0, str(script_dir.parents[2]))
from cvd_risk import PatientData
from cvd_risk.models import Prevent
import cvd_risk.utils as utils

def createPrevent(m):
	model = Prevent(m)
	model._debug = True
	model.decimal_digits = 16 # to compare results with preventr we need standard rounding
	return model



# install R tools from here:
# https://cran.r-project.org/bin/windows/Rtools/rtools45/files/rtools45-6768-6492.exe
rhome = utils.getLatestFolder(r'C:\Program Files\R', 'R-*')
if not rhome:
	raise RuntimeError('R tools missing, install ')
os.environ['R_HOME'] = rhome
os.environ['PATH'] = os.environ['R_HOME'] + r'\bin\x64;' + os.environ['PATH']

try:
    from rpy2.robjects.packages import importr
finally:
    # This second call forces rpy2 to successfully complete the initialization
    from rpy2.robjects.packages import importr
import rpy2.robjects as robjects

try:
	preventr = importr('preventr')
	install = False
except:
	install = True

if install:
	print('--- Installing preventr ---')

	# if install fails delete folder
	# C:\Users\<user>\AppData\Local\R\win-library\4.6\preventr

	binary_path = 'https://cran.r-project.org/bin/windows/contrib/4.6/preventr_0.12.2.zip'
	binary_path = r"C:\INSTALL\tracker\preventr_0.12.2.zip"

	robjects.r('options(install.packages.compile.from.source = "never")')
	utils = importr('utils')
	utils.install_packages(binary_path, repos=robjects.r('NULL'), type='win.binary')
	preventr = importr('preventr')
	print('--- Binary Installed Successfully! ---')

print('preventr version:', preventr.__version__)



def clearRes():
	global res, dbgA2, dbgB2, dbgC2
	res = collections.defaultdict(lambda: collections.defaultdict(list))
	dbgA2 = collections.defaultdict(lambda: [])
	dbgB2 = collections.defaultdict(lambda: [])
	dbgC2 = collections.defaultdict(lambda: [])

models = [
	"total_cvd_10", "ascvd_10", "hf_10", "chd_10", "stroke_10",
	"total_cvd_30", "ascvd_30", "hf_30", "chd_30", "stroke_30"
]

def calculate_row_risk(row):
	#print('='*99)
	#print(row)
	#print(help(preventr.estimate_risk))
	#xxxxxx
	if 'us_zip_code' in row and pd.notna(row['us_zip_code']):
		usz = str(int(row['us_zip_code']))
		while len(usz) < 5:
			usz = '0' + usz
	else:
		usz = robjects.NULL
	r = preventr.estimate_risk(
		age = row['age'],
		sex = row['sex'],
		sbp = row['systolic_bp'],
		bp_tx = row['antihypertensive'],
		total_c = row['total_cholesterol'] / 0.02586,
		hdl_c = row['hdl_cholesterol'] / 0.02586,
		statin = row['statin_use'],
		dm = row['diabetes'],
		egfr = row['egfr'],
		bmi = row['bmi'],
		zip = usz,
		hba1c = row['hba1c_percent'] if 'hba1c_percent' in row and pd.notna(row['hba1c_percent']) else robjects.NA_Logical,
		uacr = row['acr'] if 'acr' in row and pd.notna(row['acr']) else robjects.NA_Logical,
		smoking = row['smoking'],
		quiet = False,
	)
	r10 = r.rx2('risk_est_10yr')
	r30 = r.rx2('risk_est_30yr')
	d = {
		'total_cvd_10': r10.rx2('total_cvd'),
		'hf_10': r10.rx2('heart_failure'),
		'ascvd_10': r10.rx2('ascvd'),
		'chd_10': r10.rx2('chd'),
		'stroke_10': r10.rx2('stroke'),

		'total_cvd_30': r30.rx2('total_cvd'),
		'hf_30': r30.rx2('heart_failure'),
		'ascvd_30': r30.rx2('ascvd'),
		'chd_30': r30.rx2('chd'),
		'stroke_30': r30.rx2('stroke'),
	}
	for m in models:
		v = utils.round_half_up(d[m][0] * 100, 16)
		res[m][row['age']].append(('r', v))

	data = PatientData(
		age = row['age'],
		sex = row['sex'],
		systolic_bp = row['systolic_bp'],
		antihypertensive = row['antihypertensive'],
		total_cholesterol = row['total_cholesterol'],
		hdl_cholesterol = row['hdl_cholesterol'],
		statin_use = row['statin_use'],
		diabetes = row['diabetes'],
		egfr = row['egfr'],
		bmi = row['bmi'],
		us_zip_code = row['us_zip_code'] if 'us_zip_code' in row else None,
		hba1c_percent = row['hba1c_percent'] if 'hba1c_percent' in row and pd.notna(row['hba1c_percent']) else None,
		acr = row['acr'] if 'acr' in row and pd.notna(row['acr']) else None,
		smoking = row['smoking'],
	)
	for m in models:
		model = createPrevent(m)
		results = model.calculate(data)
		dbgA2[m].append(results.calculation_metadata['dbgA'])
		dbgB2[m].append(results.calculation_metadata['dbgB'])
		dbgC2[m].append(results.calculation_metadata['dbgC'])
		res[m][row['age']].append(('s', results.risk_score))



# regressions tests
d = {
   'age': [55, 58, 45, 49, 39],
   'sex': ['male', 'female', 'male', 'male', 'female'],
   'systolic_bp': [140.0, 130.0, 150.0, 105, 130],
   'total_cholesterol': [6.0, 5.5, 7.0, 5.6, 8],
   'hdl_cholesterol': [1.2, 1.5, 1.0, 1.1, .7],
   'smoking': [True, False, True, True, False],
   'statin_use': [True, False, True, True, True],
   'diabetes': [True, False, True, True, True],
   'antihypertensive': [True, False, True, True, True],
   'egfr': [98, 140, 60, 15, 77], # 15 <= egfr <= 140
   'bmi': [22, 23, 19, 26, 25],
   'sdi': [None, 1, None, None, 4],
   'us_zip_code': [None, None, '90813', '90808', 90814], # 1=90808, 2=90814, 3=90813
   'acr': [None, 9, None, 10, None],
   'hba1c_percent': [5.3, None, None, 5.9, None],
}

def add_testcase(d2):
	for k, v in d2.items():
		d[k] += v

# this test case causes rpy2 to return -1e-30 instead of a low risk value smaller 0.05
d2 = {'acr': [None],
 'age': [32.9],
 'antihypertensive': [False],
 'bmi': [22.7],
 'diabetes': [True],
 'egfr': [115],
 'hba1c_percent': [4.6],
 'hdl_cholesterol': [2.56],
 'sex': ['female'],
 'smoking': [False],
 'statin_use': [True],
 'systolic_bp': [94],
 'total_cholesterol': [4.69],
 'sdi': [None],
 'us_zip_code': [None]}
add_testcase(d2)


# do regression tests

clearRes()
df = pd.DataFrame(d)
# import numpy as np

# converted_zip = pd.to_numeric(df["us_zip_code"], errors="coerce").astype("Int32")
# df["sdi"] = np.where(df["sdi"].isna(), converted_zip, df["sdi"])
# print(converted_zip)
# print(df["sdi"])
# sdf

# results = df
# results["sdi"] = None
# mask = results["sdi"].isna() & results["us_zip_code"].notna()
# print(mask)
# results.loc[mask, "us_zip_code"] = results.loc[mask, "us_zip_code"].astype(str)
# print(results["us_zip_code"].astype(int))
# sdf


df.apply(calculate_row_risk, axis=1)

for m in models:
	model = createPrevent(m)
	results = model.calculate_batch(df)
	results['dbgA2'] = dbgA2[m]
	results['dbgB2'] = dbgB2[m]
	results['dbgC2'] = dbgC2[m]
	print('*'*99)
	print(m)
	print(results[['model', 'age', 'sdi', 'zip_converted', 'dbgA', 'dbgA2', 'dbgB', 'dbgB2', 'dbgC', 'dbgC2']])
	for age, risk_score in zip(results['age'], results['risk_score']):
		res[m][age].append(('b', risk_score))

nErrors = 0
for m in res.values():
	for n in m.values():
		for t, rr in n:
			if t == 'r':
				break
		ok = 1
		for t, r in n:
			if utils.round_half_up(r, 1) != utils.round_half_up(rr, 1):
				ok = 0
		if not ok:
			n.append('NOK')
			nErrors += 1

print()
print('d = ', end='')
pprint.pprint(d)
print()
pprint.pprint(res)
print(nErrors, 'nErrors')


# do random tests

for i in range(999999):

	print(i, '#'*99)
	clearRes()

	d = {
	   'age': [round(random.uniform(30, 59), 1)],
	   'sex': random.choices(['male', 'female']),
	   'systolic_bp': [random.randrange(90, 180)],
	   'total_cholesterol': [round(random.uniform(3.37, 8.27), 2)],
	   'hdl_cholesterol': [round(random.uniform(0.53, 2.58), 2)],
	   'smoking': random.choices([True, False]),
	   'statin_use': random.choices([True, False]),
	   'diabetes': random.choices([True, False]),
	   'antihypertensive': random.choices([True, False]),
	   'egfr': [random.randrange(15, 140)],
	   'bmi': [round(random.uniform(18.5, 39.9), 1)],
	   'sdi': [None],
	   'us_zip_code': [random.randrange(500, 99999) if random.randrange(2) else None],
	   'acr': [round(random.uniform(0.1, 60), 1) if random.randrange(2) else None],
	   'hba1c_percent': [round(random.uniform(4.5, 15), 1) if random.randrange(2) else None],
	}
	#d['us_zip_code'] = random.choices(prevent_data.SDI_GROUP_nan)
	df = pd.DataFrame(d)
	df.apply(calculate_row_risk, axis=1)

	for m in models:
		model = createPrevent(m)
		results = model.calculate_batch(df)
		results['dbgA2'] = dbgA2[m]
		results['dbgB2'] = dbgB2[m]
		results['dbgC2'] = dbgC2[m]
		# print('*'*99)
		# print(m)
		# print(results[['model', 'age', 'dbgA', 'dbgA2', 'dbgB', 'dbgB2', 'dbgC', 'dbgC2']])
		for age, risk_score in zip(results['age'], results['risk_score']):
			res[m][age].append(('b', risk_score))

	nErrors = 0
	for m in res.values():
		for n in m.values():
			for t, rr in n:
				if t == 'r':
					break
			for t, r in list(n):
				if utils.round_half_up(r, 1) != utils.round_half_up(rr, 1):
					if rr != -1e-30 or r > 0.05: #XXX rpy2 is returning -1e-30 for low risk values (< 0.05) ?!?!
						n.append('NOK')
						nErrors += 1
	print()
	print('d = ', end='')
	pprint.pprint(d)
	print()
	if nErrors:
		pprint.pprint(res)
		print()
		print(nErrors, 'nErrors')
		input('press enter to continue ...')
