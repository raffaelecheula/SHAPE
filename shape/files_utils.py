################################################################################
# Raffaele Cheula, LCCP, Politecnico di Milano, cheula.raffaele@gmail.com
################################################################################

################################################################################
# PRINT VECTOR TO FILE
################################################################################

def print_vector_to_file(filename, vector_name, vector, vector_type = 'list'):

    if vector_type == 'list':
        parenthesis = ['[', ']']
    elif vector_type == 'tuple':
        parenthesis = ['(', ')']

    f = open(filename, 'a+')

    f.write("\n'"+vector_name+"'\n\n")
    f.write(vector_name+' = '+parenthesis[0])
    for i in range(len(vector)):
        if i != 0:
            f.write(' '*(len(vector_name)+4))
        f.write('{:.8f}'.format(vector[i]))
        if i != len(vector)-1:
            f.write(',\n')
        else:
            f.write(parenthesis[1]+'\n\n')

    f.close()

################################################################################
# END
################################################################################
