import utils as ut
import sys

filename = sys.argv[1]
col_num = int(sys.argv[2])

ufo_list = ut.file_reader(filename, print_header=True, delimiter=',')
states = ut.col_values(ufo_list, col_num)
state_counts = ut.value_count(ufo_list, col_num, states)

print("\nState : count")
for i in state_counts:
    print(i[0], " : ", i[1])
    
