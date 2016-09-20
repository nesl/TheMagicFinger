""" This script copies data from the waiting room folder (/data/waiting_room)
To their target location based on the given action and user name.
"""

import os
import sys

if (len(sys.argv) != 3):
	sys.exit('Usage: %s action-name user-name ' %(sys.argv[0]))

action_name = sys.argv[1]
user_name = sys.argv[2]

waiting_room_files = [x for x in os.listdir('data/waiting_room') if x.startswith("wear_")]
print ("Found %d files in waiting room.\n" %(len(waiting_room_files)))

# Compute the unique experiments
waiting_room_files_prefix = set(map(lambda x: x[5:20], waiting_room_files))

target_dir = format("data/signs/%s/%s" %(action_name, user_name))
if not os.path.exists(target_dir):
	os.makedirs(target_dir)

previous_exp = [x for x in os.listdir(target_dir) ]
new_exp_idx = len(previous_exp) + 1
print (previous_exp)

# Copy each experiment to its target folder
for exp_prefix in waiting_room_files_prefix:
	exp_files = [x for x in waiting_room_files if x.find(exp_prefix) != -1]
	if (len(exp_files) == 5):
		# Create experiment target_dir
		new_exp_path = format("%s/%04d" %(target_dir, new_exp_idx))
		os.makedirs(new_exp_path)
		# Copy files to experiment target dir
		for file_name in exp_files:
			os.rename(format("data/waiting_room/%s" %file_name), format("%s/%s" %(new_exp_path, file_name)))
			print("%s ---> %s." %(format("data/waiting_room/%s" %file_name), format("%s/%s" %(new_exp_path, file_name))))
		print("")
		new_exp_idx += 1
	else:
		print("!!! Skipping incomplete experiment.\n")