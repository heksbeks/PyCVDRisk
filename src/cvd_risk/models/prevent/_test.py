import collections, importlib.util, math, os, pprint, subprocess, sys, random
from pathlib import Path

packages = ['pandas', 'numpy', 'rpy2']
missing = [pkg for pkg in packages if importlib.util.find_spec(pkg.replace('-', '_')) is None]
if missing:
	print(f'Installing missing packages {missing}...')
	subprocess.check_call([sys.executable, '-m', 'pip', 'install'] + missing)
	print('Finished installing packages!\n')

import pandas as pd

script_dir = Path(__file__).parent.resolve()
sys.path.insert(0, str(script_dir.parents[2]))
from cvd_risk import PatientData
from cvd_risk.models import prevent
import prevent_data

def getLatestFolder(root, pattern):
	dirs = [f for f in Path(root).glob(pattern) if f.is_dir()]
	if not dirs:
		return None
	return str(max(dirs, key = lambda f: f.stat().st_ctime))

# install R tools from here:
# https://cran.r-project.org/bin/windows/Rtools/rtools45/files/rtools45-6768-6492.exe
rhome = getLatestFolder(r'C:\Program Files\R', 'R-*')
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

	binary_path = 'https://cran.r-project.org/bin/windows/contrib/4.6/preventr_0.12.0.zip'
	binary_path = r"C:\INSTALL\tracker\preventr_0.12.1.zip"

	robjects.r('options(install.packages.compile.from.source = "never")')
	utils = importr('utils')
	utils.install_packages(binary_path, repos=robjects.r('NULL'), type='win.binary')
	preventr = importr('preventr')
	print('--- Binary Installed Successfully! ---')





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

def round_half_up(x, ndigits=None):
	if not math.isfinite(x):
		return x
	d = 0.5 if x > 0 else -0.5
	if ndigits is None:
		return int(math.floor(x + d))
	m = 10**ndigits
	return math.floor(x*m + d)/m

def calculate_row_risk(row):
	#print('='*99)
	#print(row)
	#print(help(preventr.estimate_risk))
	#xxxxxx
	if 'us_zip_code' in row and row['us_zip_code'] and row['us_zip_code'] > 0:
		usz = str(int(row['us_zip_code']))
		while len(usz) < 5:
			usz = '0' + usz
	else:
		usz = robjects.NULL
	r = preventr.estimate_risk(
		age = row['age'],
		sex = row['sex'],
		sbp = row['systolic_bp'],
		bp_tx = bool(row['antihypertensive']),
		total_c = row['total_cholesterol'] / 0.02586,
		hdl_c = row['hdl_cholesterol'] / 0.02586,
		statin = bool(row['statin_use']),
		dm = bool(row['diabetes']),
		egfr = row['egfr'],
		bmi = row['bmi'],
		zip = usz,
		hba1c = row['hba1c_percent'] if 'hba1c_percent' in row and pd.notna(row['hba1c_percent']) else robjects.NA_Logical,
		uacr = float(row['acr']) if 'acr' in row and pd.notna(row['acr']) else robjects.NA_Logical,
		smoking = int(row['smoking']),
		quiet = False,
	)
	r10 = r.rx2('risk_est_10yr')
	r30 = r.rx2('risk_est_30yr')
	n = 30
	d = {
		'total_cvd_10': round_half_up(r10.rx2('total_cvd')[0] * 100, n),
		'hf_10': round_half_up(r10.rx2('heart_failure')[0] * 100, n),
		'ascvd_10': round_half_up(r10.rx2('ascvd')[0] * 100, n),
		'chd_10': round_half_up(r10.rx2('chd')[0] * 100, n),
		'stroke_10': round_half_up(r10.rx2('stroke')[0] * 100, n),

		'total_cvd_30': round_half_up(r30.rx2('total_cvd')[0] * 100, n),
		'hf_30': round_half_up(r30.rx2('heart_failure')[0] * 100, n),
		'ascvd_30': round_half_up(r30.rx2('ascvd')[0] * 100, n),
		'chd_30': round_half_up(r30.rx2('chd')[0] * 100, n),
		'stroke_30': round_half_up(r30.rx2('stroke')[0] * 100, n),
	}
	for m in models:
		res[m][row['age']].append(('r', d[m]))
	#print(d)

	data = PatientData(
		age = row['age'],
		sex = row['sex'],
		systolic_bp = row['systolic_bp'],
		antihypertensive = bool(row['antihypertensive']),
		total_cholesterol = row['total_cholesterol'],
		hdl_cholesterol = row['hdl_cholesterol'],
		statin_use = bool(row['statin_use']),
		diabetes = bool(row['diabetes']),
		egfr = row['egfr'],
		bmi = row['bmi'],
		us_zip_code = int(row['us_zip_code']) if 'us_zip_code' in row and row['us_zip_code'] and row['us_zip_code'] > 0 else None,
		hba1c_percent = row['hba1c_percent'] if 'hba1c_percent' in row and pd.notna(row['hba1c_percent']) else None,
		acr = float(row['acr']) if 'acr' in row and pd.notna(row['acr']) else None,
		smoking = int(row['smoking']),
	)
	for m in models:
		model = prevent.Prevent(m)
		model._debug = True
		model._round = lambda x: round_half_up(x, 30)
		results = model.calculate(data)
		dbgA2[m].append(results.calculation_metadata['dbgA'])
		dbgB2[m].append(results.calculation_metadata['dbgB'])
		dbgC2[m].append(results.calculation_metadata['dbgC'])
		res[m][row['age']].append(('s', results.risk_score))



# Create DataFrame with patient data
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
   #'sdi': [None, None, None, None, None],
   'us_zip_code': [None, None, 90813, 90815, None], # 1=90815, 2=90814, 3=90813
   'acr': [None, 9, None, 10, None],
   'hba1c_percent': [5.3, None, None, 5.9, None],
}
# for k, v in d.items():
	# d[k] = [v[3]]

d2 = {'acr': [None],
 'age': [39.0],
 'antihypertensive': [False],
 'bmi': [22.1],
 'diabetes': [False],
 'egfr': [139],
 'hba1c_percent': [4.6],
 'hdl_cholesterol': [2.57],
 'sex': ['female'],
 'smoking': [False],
 'statin_use': [True],
 'systolic_bp': [180],
 'total_cholesterol': [3.37],
 'us_zip_code': [None]}

df = pd.DataFrame(d)


# Calculate risks

clearRes()
df.apply(calculate_row_risk, axis=1)

for m in models:
	model = prevent.Prevent(m)
	model._debug = True
	results = model.calculate_batch(df)
	results['dbgA2'] = dbgA2[m]
	results['dbgB2'] = dbgB2[m]
	results['dbgC2'] = dbgC2[m]
	print('*'*99)
	print(m)
	print(results[['model', 'age', 'dbgA', 'dbgA2', 'dbgB', 'dbgB2', 'dbgC', 'dbgC2']])
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
			if round_half_up(r, 1) != round_half_up(rr, 1):
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
	   #'sdi': [None, None, None, None, None],
	   'us_zip_code': [random.randrange(500, 99999) if random.randrange(2) else None],
	   'acr': [round(random.uniform(0.1, 60), 1) if random.randrange(2) else None],
	   'hba1c_percent': [round(random.uniform(4.5, 15), 1) if random.randrange(2) else None],
	}
	#d['us_zip_code'] = random.choices(prevent_data.SDI_GROUP_nan)
	df = pd.DataFrame(d)
	df.apply(calculate_row_risk, axis=1)

	for m in models:
		model = prevent.Prevent(m)
		model._debug = True
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
				if round_half_up(r, 1) != round_half_up(rr, 1):
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
