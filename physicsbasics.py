
def standard_form(array,sf): 
    for i in range(len(array)):
        number = array[i]
        n = 0
        while number>10:
            number = number/10
            n = n+1
        while number<1: 
            number = number*10
            n = n-1
        array[i] = (round(number,sf))*10**n
    return array


