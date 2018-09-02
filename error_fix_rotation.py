import numpy as np
import itertools
import matplotlib.pyplot as plt


def bin_to_num(binary):
    num = np.sum([2**-(n+1) for n,i in enumerate(reversed(binary)) if i=='1'])
    return num

def max_sign_bit(binary):
    return len(binary)-list(reversed(binary)).index('1')

def min_lamb(k):
    return 2**-k
def get_error(lamb,msb,accuracy,k,choose):
    if choose == 'max':
        er_lamb = 1/2**(k-msb+accuracy+1)# error of lambda = 0.5 bin width
        er_lamb[msb<=accuracy] = 0
        max_error = np.max(np.abs((np.arcsin(min_lamb(k) / (lamb - er_lamb)) -
                                   np.arcsin(min_lamb(k) / lamb))))
        print(lamb,er_lamb,msb)
        return max_error
    elif choose == 'maxrel':
        er_lamb = 1/2**(k-msb+accuracy+1)# error of lambda = 0.5 bin width
        er_lamb[msb<=accuracy] = 0
        maxrel_error = np.max(np.abs((np.arcsin(min_lamb(k) / (lamb - er_lamb)) -
                                   np.arcsin(min_lamb(k) / lamb)) /
                                  np.arcsin(min_lamb(k) / (lamb - er_lamb))))
        #print(lamb,er_lamb,msb)
        return maxrel_error
    
    else:
        raise ValueError('{} not yet implemented'.format(choose))
def get_error_2(lamb,msb,accuracy,k,choose):
    if choose == 'max':
        er_lamb = 1/2**(k-msb+accuracy+1)# error of lambda = 0.5 bin width
        er_lamb[msb<=accuracy] = 0
        max_error = np.max(np.abs(((min_lamb(k) / (lamb - er_lamb)) -
                                   (min_lamb(k) / lamb))))
        #print(lamb,er_lamb,msb)
        return max_error
    elif choose == 'maxrel':
        er_lamb = 1/2**(k-msb+accuracy+1)# error of lambda = 0.5 bin width
        er_lamb[msb<=accuracy] = 0
        maxrel_error = np.max(np.abs(((min_lamb(k) / (lamb - er_lamb)) -
                                   (min_lamb(k) / lamb)) /
                                  (min_lamb(k) / (lamb - er_lamb))))
        #print(lamb,er_lamb,msb,accuracy)
        return maxrel_error
    
    else:
        raise ValueError('{} not yet implemented'.format(choose))
def get_est_lamb(pattern,msb,n):
    '''Estimate the bin mid point and return the float value'''
    if msb-n > 0:
        pattern[msb-n-1]='1'
        return bin_to_num(pattern)
    else:
        return bin_to_num(pattern)
    
def error_analysis(k,n):
    '''Calculate error of arcsin rotation using k bits fixed point numbers and n bit accuracy'''
    lambda_array = []
    msb_array = []
    n_ = n #copy of n
    for msb in range(k-1,-1,-1):
        vec = ['0']*k
        vec[msb] = '1'
        if msb <= n_-1:
            n -= 1
            #print(msb,n,n_)
        for pattern in itertools.product('10',repeat=n):
            vec[msb-n:msb] = pattern
            e_l = get_est_lamb(vec.copy(),msb,n)
            l = bin_to_num(vec)
            #print(vec,e_l,l,get_error_2(np.array([e_l]),np.array([msb]),n_,k,'maxrel'))
            lambda_array.append(get_est_lamb(vec.copy(), msb, n))
            msb_array.append(msb)
    #print("finished here")
    return (get_error_2(np.array(lambda_array),np.array(msb_array),n_,k,'maxrel'))
            
#error_analysis(6,4)
def newton_iteration(lam,acu):
    x0 = 2**-np.ceil(np.log2(lam))
    #print(x0)
    def get_storage(lam,acu):
        _ = 0
        lam_orig = lam
        bin_ = bin(int(lam))
        aprox_lam = int(bin_,2)
        while abs(aprox_lam-lam)/lam > acu:
            _ += 1
            lam *= 10
            bin_ = bin(int(lam))
            aprox_lam = int(bin_,2)
            if _ >= 10:
                raise RuntimeError("Cannot store inverse eigenvalue")
        #print(lam_orig,lam,aprox_lam,_)
        return _ #np.floor(np.log2(lam))+1
            

    def error(lam,x,thres = 0.05):
        if np.abs(1/lam-x)/(1/lam) > 0.05:
            #print(x,np.abs(1/lam-x)/(1/lam))
            return True
        else:
            return False
    _ = 0
    x = x0
    while error(lam,x,acu):
        if _ > 50:
            raise RuntimeError('Newton iteration failed')
        _ += 1
        x = - lam * x**2 + 2 * x
    #print("Number of iterations:",_)
    
    return _,get_storage(x,acu)
"""    
for k in range(1,12):
    res = []
    print('\n\n\n')
    for n in range(3,k):
        res.append(error_analysis(k,n))
        print('\n\n\n')
    plt.plot(range(3,k),res,label=str(k))
    #print(k,res)
plt.legend(loc='best')
plt.show()
"""
max_it = []
max_bit = []
acu = np.sqrt(0.05**2 / 2)
for k in range(4,12):
    res_it = []
    res_bit = []
    for pattern in itertools.product('10',repeat=k):
        lam = bin_to_num(pattern)
        if lam == 0:
            continue
        iterations,bits = newton_iteration(lam,acu)
        res_it.append(iterations)
        res_bit.append(bits+len(pattern)+1)
    max_it.append(np.max(res_it))
    max_bit.append(np.max(res_bit))
    print("Minimum iterations to get accuracy of {} for k = {} : {}".format(acu,k,np.max(res_it)))
    print("Storage needed: {} bits".format(np.max(res_bit)))
