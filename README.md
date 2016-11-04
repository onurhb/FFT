# `FFT`
`Simple and efficent FFT implementation`

**`When debugging butterfly algorithm it should output something like this for 1hz sine, N = 8`**

`k : indicates where butterfly starts (upper)
 b_size : size of the butterfly
 b_count : butterfly offset`
 
```
Stage : 1
Butterfly -> k: 0 b_size: 1 b_count: 4
 ------------------------------
k0 = k0 + k1 * W^0_2
k1 = k1 - k0 * W^0_2
------------------------------
Butterfly -> k: 2 b_size: 1 b_count: 4
 ------------------------------
k2 = k2 + k3 * W^0_2
k3 = k3 - k2 * W^0_2
------------------------------
Butterfly -> k: 4 b_size: 1 b_count: 4
 ------------------------------
k4 = k4 + k5 * W^0_2
k5 = k5 - k4 * W^0_2
------------------------------
Butterfly -> k: 6 b_size: 1 b_count: 4
 ------------------------------
k6 = k6 + k7 * W^0_2
k7 = k7 - k6 * W^0_2
------------------------------
Stage : 2
Butterfly -> k: 0 b_size: 2 b_count: 2
 ------------------------------
k0 = k0 + k2 * W^0_4
k2 = k2 - k0 * W^0_4
------------------------------
Butterfly -> k: 1 b_size: 2 b_count: 2
 ------------------------------
k1 = k1 + k3 * W^1_4
k3 = k3 - k1 * W^1_4
------------------------------
Butterfly -> k: 4 b_size: 2 b_count: 2
 ------------------------------
k4 = k4 + k6 * W^0_4
k6 = k6 - k4 * W^0_4
------------------------------
Butterfly -> k: 5 b_size: 2 b_count: 2
 ------------------------------
k5 = k5 + k7 * W^1_4
k7 = k7 - k5 * W^1_4
------------------------------
Stage : 3
Butterfly -> k: 0 b_size: 4 b_count: 1
 ------------------------------
k0 = k0 + k4 * W^0_8
k4 = k4 - k0 * W^0_8
------------------------------
Butterfly -> k: 1 b_size: 4 b_count: 1
 ------------------------------
k1 = k1 + k5 * W^1_8
k5 = k5 - k1 * W^1_8
------------------------------
Butterfly -> k: 2 b_size: 4 b_count: 1
 ------------------------------
k2 = k2 + k6 * W^2_8
k6 = k6 - k2 * W^2_8
------------------------------
Butterfly -> k: 3 b_size: 4 b_count: 1
 ------------------------------
k3 = k3 + k7 * W^3_8
k7 = k7 - k3 * W^3_8
------------------------------
0 Hz : 0 dbs
1 Hz : 1 dbs
2 Hz : 0 dbs
3 Hz : 0 dbs
Time: 0.008000
```


![Butterfly diagram](https://github.com/onurmhb/FFT/blob/master/butterfly.png "Butterfly")
