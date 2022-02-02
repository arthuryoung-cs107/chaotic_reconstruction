clear
close all

sym2 = [1.000000 2.000000 3.000000 4.000000; 2.000000 5.000000 6.000000 7.000000; 3.000000 6.000000 8.000000 9.000000; 4.000000 7.000000 9.000000 10.000000]

m1 = [1.000000 2.000000 3.000000 4.000000; 5.000000 6.000000 7.000000 8.000000; 9.000000 10.000000 11.000000 12.000000; 13.000000 14.000000 15.000000 16.000000; 17.000000 18.000000 19.000000 20.000000];

v1 =  (1:4)'
v3 = sym2*v1 + v1
v4 = v1 - sym2*v1
vT_A_v = v1'*sym2*v1

m1
sym3 = m1'*m1

sym4 =  [1.357711 0.621911 0.656543 0.379332; 0.621911 0.909706 0.507339 0.987917; 0.656543 0.507339 0.589831 0.635260; 0.379332 0.987917 0.635260 1.444418]
v6 = [0.716341; 0.362693; 0.475673; 0.272259]
v7 = [0.105670; 0.004159; 0.944203; 0.335002]

check = v6'*(sym4*v6 + v7)
v7 = sym4*v6 + v7
