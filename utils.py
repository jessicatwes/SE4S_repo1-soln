def file_reader(filename, print_header=True, delimiter=','):
    '''
    This function reads in the entire csv file into a list of lists.
    '''
    ufo_list = []
    with open(filename) as ufo_data:
        header = ufo_data.readline().strip().split(delimiter)

        for line in ufo_data:
            line = line.strip().split(delimiter)
            ufo_list.append(line)

    if print_header == True:
        print("File header: ", ",".join(header))

    return ufo_list


def col_values(data, col_num):
    '''
    This function takes the list of lists from file_reader() and returns a 
    list of all the unique entries for a given column.
    '''
    val_set = set()
    for row in data:
        val_set.add(row[col_num])
    
    values = list(val_set)
    return list(values)


def value_count(data, col_num, values):
    '''
    This function takes the data from file_reader() and the list of values
    from col_values() and then returns a list of len-2 lists with the
    count of each value in the original data. 
    '''
    count_list = [[i, 0] for i in sorted(values)]
    for row in data:
        row_val = row[col_num]

        for i, val_count in enumerate(count_list):

            if val_count[0] == row_val:
                count_list[i][1] += 1
                break
            else:
                continue

    return count_list
