""" Load data from CSV files
"""

import os
import numpy as np

def load_exp(sign, user, exp_id):
	exp_path = format("data/signs/%s/%s/%s" %(sign, user, exp_id))
	exp_files = [x for x in os.listdir(exp_path) if x.endswith(".csv")]
	acc_data = []
	gyro_data = []
	mag_data = []
	grav_data = []
	rot_data = []

	for file_name in exp_files:
		if (file_name.endswith("acc.csv")):
			acc_data = np.loadtxt(format("%s/%s" %(exp_path, file_name)), delimiter=',')
		elif (file_name.endswith("gyro.csv")):
			gyro_data = np.loadtxt(format("%s/%s" %(exp_path, file_name)), delimiter=',')
		elif (file_name.endswith("mag.csv")):
			mag_data = np.loadtxt(format("%s/%s" %(exp_path, file_name)), delimiter=',')
		elif (file_name.endswith("rot.csv")):
			rot_data = np.loadtxt(format("%s/%s" %(exp_path, file_name)), delimiter=',')
		elif (file_name.endswith("grav.csv")):
			grav_data = np.loadtxt(format("%s/%s" %(exp_path, file_name)), delimiter=',')

	exp_data = {
	'id': exp_id,
	'user': user,
	'sign': sign,
	'acc': acc_data,
	'gyro': gyro_data,
	'mag': mag_data,
	'grav': grav_data,
	'rot': rot_data}
	return exp_data

def load_all_experiments():
    experiments_list = []
    signs = [x for x in os.listdir("data/signs") if os.path.isdir(os.path.join("data/signs", x))]
    for sign_name in signs:
        sign_path = format("data/signs/%s" %sign_name)
        users_list = [x for x in os.listdir(sign_path) if os.path.isdir(os.path.join(sign_path, x))]
        for user_name in users_list:
            user_path = os.path.join(sign_path, user_name)
            exp_list = [x for x in os.listdir(user_path) if os.path.isdir(os.path.join(user_path, x))]
            for exp_id in exp_list:
                exp_data = load_exp(sign_name, user_name, exp_id)
                experiments_list.append(exp_data)
    return experiments_list
