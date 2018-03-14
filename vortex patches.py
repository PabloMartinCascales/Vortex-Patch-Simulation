# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 20:33:25 2018

@author: morot
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Jan 14 20:28:12 2018

@author: morot
"""
import time
import numpy as np
import math
A = 2.
B = 1.
#PERIMETRO = 2*PI*math.sqrt(((A**2)+(B**2))/2.)
#Para A y B dados
PERIMETRO = 9.11378427674
PERIMETRO_2 = 28.16976286
TOTAL_ITER = 500
DT = 0.005

PI = math.pi


def simpson_compuesto(vec, elem, a, b, n):
    #	calculamos primero h con el número n de subintervalos, la parametrizacion ya la tenemos en el vector
    # para un vector de dim 100: b=99, a=0.
    dim = 2*n
    #h = ((8*(b+1)/(dim))-(a/(dim)))/(dim)
    h = ((8.*(b+1)/(dim))-(a/(dim)))/(dim)
    
    impar = 0.
    par = 0.
    
    for idx_1 in range(1, n):
        
        impar = impar + f(vec, elem, a+(2*idx_1-1), dim)
        #if (2*idx_1-1) == b:
           #print('wachi')
            
    for idx_2 in range(1, n):
        par = par + f(vec, elem, a+(2*idx_2), dim)
    # esto funciona
    return (1/(2*math.pi))*((h/3.)*(f(vec, elem, a, dim) + f(vec, elem, b, dim) + 2*par + 4*impar))

def f(vect , e, x, dim):
    # s denota la posicion en el vector parametrizado vector y x es la variable muda de la integración,    
    
    dif = e-vect[x]
    if (list(dif) == [0, 0]):
        return 0.
    else:
        return ((dif)/((np.linalg.norm(dif))**2))*np.dot(dif, f_aux(vect , x, dim))

#print simpson_compuesto(-1., 4., NUM_SUB)
def f_aux(vect, x, dim):
    
    
    if x+1 == dim:
        dif = vect[0+1]-vect[x]
        return (dif)/(np.linalg.norm(dif))
    elif x+1 < dim:
        dif = vect[x+1]-vect[x]
        return (dif)/(np.linalg.norm(dif))
    
    
    
    
    

def runge_kutta_4(vect, x_0, dim):
    # Aquí vamos a tratar una edo de la forma x'(t) = g(t,x), tenemos que dar la condición incial, el valor de vector[s] en un paso temporal dado
    # Solo hacemos el runge kutta una vez hata que actualicemos todos los puntos del vector
    n_inter=int(dim/2)
    F_1 = DT*simpson_compuesto(vect, x_0, 0, dim-1, n_inter)
    
    for dummy_1 in range(2):    
        F_2 = DT*simpson_compuesto(vect, x_0+(F_1/2.), 0, dim-1, n_inter)
        F_3 = DT*simpson_compuesto(vect, x_0+(F_2/2.), 0, dim-1, n_inter)
        F_4 = DT*simpson_compuesto(vect, x_0+F_3, 0, dim-1, n_inter)
        
        x_1 = x_0 + F_1/6. + F_2/3. + F_3/3. + F_4/6.
        x_0 = x_1

    return x_1


def area_vector(vec, a, b, n, function):
    #	calculamos primero h con el número n de subintervalos, la parametrizacion ya la tenemos en el vector
    # para un vector de dim 100: b=99, a=0.
    #h = ((4*(b+1)/(2*n))-(a/(2*n)))/(2*n)
    h = (((b+1)/(2*n))-(a/(2*n)))/(2*n)
    impar = 0.
    par = 0.
    
    for idx_1 in range(1, n):
        
        impar = impar + function(vec, a+(2*idx_1-1))
        #if (2*idx_1-1) == b:
           #print('wachi')
            
    for idx_2 in range(1, n):
        par = par + function(vec, a+(2*idx_2))
    # Para el circulo multiplicamos por 2pi
    #return (((h/3.)*(2*PI)*(function(vec, a) + function(vec, b) + 2*par + 4*impar))/2.)
    #return (((h/3.)*PERIMETRO*(function(vec, a) + function(vec, b) + 2*par + 4*impar))/2.)    
    return (((h/3.)*PERIMETRO*(function(vec, a) + function(vec, b) + 2*par + 4*impar))/2.)


def f_general(vec, x):
    s = normal_vector(tangente(vec, x))
    return np.dot(vec[x], s)
#Aqui ponemos en tangente una dimension de 500 del vector para ahoorar t computacion
def tangente(vec, i):
    if i == 0:
        return (vec[i+1]-vec[499])
    elif i == 499:
        return (vec[0]-vec[i-1])
    else:
        return (vec[i+1]-vec[i-1])


def normal_vector(tan):
    z = np.array([0,0,1])
    v = np.append(tan,[0])
    
    normal = np.cross(v,z)
    normal = normal/(np.linalg.norm(normal))
    
    return (np.delete(normal,2))


def f_circle(vec, x):
    return np.linalg.norm(vec[x])


def check_vector(vector, tol):
    
    dim = len(vector)
    
    aux_vector = vector
    i = 0
    while i+2 < dim:   
        dif_1 = aux_vector[i+1]-aux_vector[i]
        dif_2 = aux_vector[i+1]-aux_vector[i+1]
        dif_3 = aux_vector[i+2]-aux_vector[i]
        dist_1 = np.linalg.norm(dif_1)
        dist_2 = np.linalg.norm(dif_2)
        
        if ((dist_1 <= tol) and (dist_2 > tol)) or (dist_1 > tol) and (dist_2 <= tol):
            
            aux_vector[i+1] = aux_vector[i] + dif_3/2. 
        
        elif (dist_1 > tol) and (dist_2 > tol):
            aux_vector[i+1] = aux_vector[i] + dif_1/2.
            aux_vector[i+2] = aux_vector[i] + dif_3/2.
        
        i = i+1
    #Saliendo del bucle tenemos un i+2 = dim
    dif_1 = aux_vector[i+1]-aux_vector[i]
    dif_2 = aux_vector[i+1]-aux_vector[i+1]
    dif_3 = aux_vector[0]-aux_vector[i]
    dist_1 = np.linalg.norm(dif_1)
    dist_2 = np.linalg.norm(dif_2)
    
    if ((dist_1 <= tol) and (dist_2 > tol)) or (dist_1 > tol) and (dist_2 <= tol):
            
        aux_vector[i+1] = aux_vector[i] + dif_3/2.
            
        
    elif (dist_1 > tol) and (dist_2 > tol):
        aux_vector[i+1] = aux_vector[i] + dif_1/2.
        aux_vector[0] = aux_vector[i] + dif_3/2.
    
    i = i+1
    
    dif_1 = aux_vector[0]-aux_vector[i]
    dif_2 = aux_vector[1]-aux_vector[0]
    dif_3 = aux_vector[1]-aux_vector[i]
    dist_1 = np.linalg.norm(dif_1)
    dist_2 = np.linalg.norm(dif_2)
    
    if ((dist_1 <= tol) and (dist_2 > tol)) or (dist_1 > tol) and (dist_2 <= tol):
            
        aux_vector[0] = aux_vector[i] + dif_3/2.
            

    elif (dist_1 > tol) and (dist_2 > tol):
        aux_vector[0] = aux_vector[i] + dif_1/2.
        aux_vector[1] = aux_vector[i] + dif_3/2.
    
    return aux_vector

def write_file(vector,name):
    file = open(name, "w")

    for idx in range(len(vector)):
        file.write(str(vector[idx][0])+" ")
        file.write(str(vector[idx][1])+"\n")
    
    file.close()
    
def read_file(name):
    file = open(name, "r")
    vector_aux = []
    a =' '
    b =' '
   
    
    for line in file:
        j = 0
        while line[j] != " ":
            a = a + line[j]
            j = j+1
        
        for k in range(j, len(line), 1):
            b = b + line[k]

        s = float(a)
        m = float(b)
        a = ' '
        b = ' '
        
        vector_aux.append(np.array([s, m]))
    
    file.close()
    return vector_aux    
    
    
def write_area(area,name):
    file = open(name, "w")
    file.write(str(area))
    file.close()

def distance(vector):
    
    lista = []
    for i in range(len(vector)-1):    
        lista.append(np.linalg.norm(vector[i+1]-vector[i])) 
    i = i+1
    
    lista.append(np.linalg.norm(vector[0]-vector[i]))
    return lista
    
    
    

###
def create_elipse():
    global tol_elipse
    Xi = np.arange(0, 1, 0.002)
    vector_aux = []
    
    for num in Xi: 
        vector_aux.append(np.array([A*math.cos(2*PI*num),B*math.sin(2*PI*num)]))
        
    lista = []
    for i in range(len(vector_aux)-1):    
        lista.append(np.linalg.norm(vector_aux[i+1]-vector_aux[i])) 
    i = i+1
    
    lista.append(np.linalg.norm(vector_aux[0]-vector_aux[i]))
    tol_elipse = max(lista)
    
    return vector_aux



def create_circle():
    global tol_circle
    Gi = np.arange(0, 1, 0.002)
    vector_aux = []

    for num in Gi: 
        vector_aux.append(np.array([math.cos(2*PI*num),math.sin(2*PI*num)]))
    
    tol_circle = np.linalg.norm(vector_aux[0]-vector_aux[1])
    return vector_aux


def create_triangle():
    # 600 puntos
    global tol_triangle
    Gi = np.arange(0, 1, 1./600.)
    vector_aux = []

    for idx_1 in range(0, 150, 1): 
        vector_aux.append(np.array([5.-(20.)*Gi[idx_1],(20.)*Gi[idx_1]]))
        
    for idx_2 in range(150, 300, 1): 
        vector_aux.append(np.array([5.-(20.)*Gi[idx_2],10.-(20.)*Gi[idx_2]]))
        
    for idx_3 in range(300, 600, 1): 
        vector_aux.append(np.array([-15.+(20.)*Gi[idx_3],0]))
        
    lista = []
    for i in range(len(vector_aux)-1):    
        lista.append(np.linalg.norm(vector_aux[i+1]-vector_aux[i])) 
    i = i+1
    
    lista.append(np.linalg.norm(vector_aux[0]-vector_aux[i]))
    tol_triangle = max(lista)
 
    return vector_aux
    
###

"""
def run_simulation(vector, times): 
    
    #Debemos poner tol_elipse o tol_circle en funcion de que simulemos   
    max_distance = max(distance(vector))
    dim = len(vector)   
    for temp in range(1, times+1):  
        
        for idx in range(len(vector)):
            
            vector[idx] = runge_kutta_4(vector, vector[idx], dim) 
            if temp%2 == 0:
                count = 1
                while max_distance > 1.1*tol_elipse:
                    vector = check_vector(vector, tol_elipse)
                    if count%5 == 0:
                        max_distance = max(distance(vector))
                    count += 1     
    return vector
"""

def run_simulation(vector, times): 
    
    #Debemos poner tol_elipse o tol_circle en funcion de que simulemos, este run simulation escribe ficheros no devuelve un vector   
    max_distance = max(distance(vector))
    dim = len(vector)   
    for temp in range(1, times+1):  
        
        for idx in range(dim):          
            vector[idx] = runge_kutta_4(vector, vector[idx], dim) 
            
        if temp%50 == 0:
            count = 1
            while max_distance > 1.1*tol_triangle:
                vector = check_vector(vector, tol_triangle)
                if count%5 == 0:
                    max_distance = max(distance(vector))
                count += 1
               
        if temp == 100:
            write_file(vector, 'tri_100.txt')
            write_area(area_vector(vector, 0, dim-1, int(dim/2), f_general), 'area_100_tri.txt')
        if temp == 200:
            write_file(vector, 'tri_200.txt')
            write_area(area_vector(vector, 0, dim-1, int(dim/2), f_general), 'area_200_tri.txt')
        elif temp == 300:
            write_file(vector, 'tri_300.txt')
            write_area(area_vector(vector, 0, dim-1, int(dim/2), f_general), 'area_300_tri.txt') 
        elif temp == 400:
            write_file(vector, 'tri_400.txt')
            write_area(area_vector(vector, 0, dim-1, int(dim/2), f_general), 'area_400_tri.txt')
        elif temp == 500:
            write_file(vector, 'tri_500.txt')
            write_area(area_vector(vector, 0, dim-1, int(dim/2), f_general), 'area_500_tri.txt')            
        elif temp == 600:
            write_file(vector, 'tri_600.txt')
            write_area(area_vector(vector, 0, dim-1, int(dim/2), f_general), 'area_600_tri.txt')            
        elif temp == 700:
            write_file(vector, 'tri_700.txt')
            write_area(area_vector(vector, 0, dim-1, int(dim/2), f_general), 'area_700_tri.txt')            
        elif temp == 800:
            write_file(vector, 'tri_800.txt')
            write_area(area_vector(vector, 0, dim-1, int(dim/2), f_general), 'area_800_tri.txt')           
        elif temp == 900:
            write_file(vector, 'tri_900.txt')
            write_area(area_vector(vector, 0, dim-1, int(dim/2), f_general), 'area_900_tri.txt')           
        elif temp == 1000:
            write_file(vector, 'tri_1000.txt')
            write_area(area_vector(vector, 0, dim-1, int(dim/2), f_general), 'area_1000_tri.txt')
        





tiempo_inicial = time.time()
###  PROGRAMA PRINCIPAL ###



b = create_elipse()
dim = len(b)
print(area_vector(b, 0, dim-1, int(dim/2), f_general))



#write_file(run_simulation(create_elipse(), 10),'x5DT.txt')
#run_simulation(create_triangle(), 1000)


"""
a = create_triangle()


for idx in range(5, 10, 1):
    a[idx] = a[idx+5]
    
write_file(a,'triangulo_roto.txt')
max_distance = max(distance(a))
print(max_distance)
print(tol_elipse)

count = 1
while max_distance > 1.1*tol_triangle:
    a = check_vector(a, tol_triangle)
    if count%5 == 0:
        max_distance = max(distance(a))
    count += 1
        
write_file(a,'triangulo_arreglado.txt')
print(count)

"""


### FIN PROGRAMA PRINCIPAL ###
tiempo_final = time.time()
tiempo_ejecucion = tiempo_final-tiempo_inicial
print('Tiempo:'+str(tiempo_ejecucion))

