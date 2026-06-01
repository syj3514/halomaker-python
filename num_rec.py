import numpy as np
from numba import njit, prange

@njit
def assign_struct_ids(nnodes, mothers, mostmasssub):
    node_to_struct = np.full(nnodes, -1, dtype=np.int32)
    nb_halos = 0
    nb_sub = 0
    istruct = 0
    
    for i in range(nnodes):
        inode1 = i + 1
        imother1 = mothers[i]
        
        if imother1 <= 0:
            nb_halos += 1
            istruct += 1
            node_to_struct[i] = istruct
        else:
            if inode1 == mostmasssub[imother1 - 1]:
                # 상속 (이 부분 때문에 단순 cumsum만으로는 안 됨)
                node_to_struct[i] = node_to_struct[imother1 - 1]
            else:
                nb_sub += 1
                istruct += 1
                node_to_struct[i] = istruct
    return node_to_struct, nb_halos, nb_sub

@njit
def counting_argsort_8(keys):
    n = keys.size
    cnt = np.zeros(8, np.int32)
    for i in range(n):
        cnt[keys[i]] += 1

    start = np.empty(8, np.int32)
    s = 0
    for k in range(8):
        start[k] = s
        s += cnt[k]
        cnt[k] = 0  # reuse as "filled"

    order = np.empty(n, np.int32)
    for i in range(n):  # stable
        k = keys[i]
        j = start[k] + cnt[k]
        order[j] = i
        cnt[k] += 1

    return order

# For calculation of the age of the universe
#***********************************************************************
def polint(xa, ya, x):
    n = len(xa)
    if n > 10:
        raise ValueError("Maximum allowed size for xa and ya is 10")
    
    c = ya.copy()
    d = ya.copy()

    # Find the index of the closest value in xa to x
    ns = np.abs(xa - x).argmin()
    y = ya[ns]
    
    for m in range(1, n):
        den = xa[:n-m] - x - xa[m:] + x
        w = c[1:n-m+1] - d[:n-m]
        with np.errstate(divide='raise'):
            try:
                den = w / den
            except FloatingPointError:
                raise ValueError("Denominator is zero!")

        d[:n-m] = (xa[m:] - x) * den
        c[:n-m] = (xa[:n-m] - x) * den

        if 2 * ns < n - m:
            dy = c[ns]
        else:
            dy = d[ns-1]
            ns -= 1
        
        y += dy
        
    return y, dy

#***********************************************************************
def qromo(func, a, b, **kwargs):
    def midpnt(func, a, b, s, n, **kwargs):
        # If `it` was a static or saved variable in Fortran subroutine,
        # then in Python, it can be made an attribute of the function to mimic this behavior.
        if not hasattr(midpnt, "it"):
            midpnt.it = 0

        if n == 1:
            s = (b - a) * func(0.5 * (a + b), **kwargs)
            midpnt.it = 1
        else:
            tnm = midpnt.it
            del_ = (b - a) / (3.0 * tnm)  # using del_ instead of del because del is a keyword in Python
            ddel = del_ + del_
            x = a + 0.5 * del_
            sum_ = 0.0  # using sum_ instead of sum because sum is a built-in function in Python

            for j in range(midpnt.it):
                sum_ += func(x, **kwargs)
                x += ddel
                sum_ += func(x, **kwargs)
                x += del_
            
            s = (s + (b - a) * sum_ / tnm) / 3.0
            midpnt.it *= 3

        return s

    jmax = 14; k = 5; km = k - 1; jmaxp = jmax + 1; dss=0
    eps = 1e-6
    s = np.zeros(jmaxp)
    h = np.zeros(jmaxp)
    h[0] = 1.0
    for j in range(jmax):
        s[j] = midpnt(func, a, b, s[j], j+1, **kwargs)
        if j >= k:
            ss, dss = polint(h[j-km : j-km+k], s[j-km : j-km+k], 0.0)
            if abs(dss) < eps * abs(ss):
                return ss
        s[j+1] = s[j]
        h[j+1] = h[j] / 9.0
    print('> too many steps in qromo')
    return ss

#***********************************************************************
def age_temps(x,omm=0,oml=0):
    # Maybe this should be moved to `compute_halo_props`
    return 1/np.sqrt(omm*x**5 + (1.0-omm-oml)*x**4 + oml*x**2)

#***********************************************************************
def age_temps_turn_around(x,omm=0,oml=0):
    # Maybe this should be moved to `compute_halo_props`
    temps = omm*x**5 + (-omm-oml)*x**4 + oml*x**2
    if(temps > 0.0):
       return 1/np.sqrt(temps)
    else:
       return 0.0

#***********************************************************************
# subroutine `locate` is not used in the code

#***********************************************************************
def jacobi(a, nmax=500, itermax=50):
    n = a.shape[0]
    assert n == a.shape[1]
    if(n>nmax): raise BufferError(f"n(={n}) is too small (nmax={nmax})")
    b = np.diag(a).copy()
    z=np.zeros(n)

    d = np.diag(a).copy()
    v = np.eye(n)
    nbinot = 0
    for i in range(itermax):
        sm = np.sum(np.abs(np.triu(a, k=1)))
        if(sm == 0): return d, v
        tresh = 0.2*sm/(n**2) if(i<4) else 0.0

        for ip in range(n-1):
            for iq in range(ip+1,n):
                g = 100*abs(a[ip,iq])
                if (i > 4) and (abs(d[ip])+g == abs(d[ip])) and (abs(d[iq])+g == abs(d[iq])):
                    a[ip,iq] = 0.0
                elif (abs(a[ip,iq]) > tresh):
                    h = d[iq] - d[ip]
                    if (abs(h)+g == abs(h)):
                        t = a[ip,iq]/h
                    else:
                        theta = 0.5*h/a[ip,iq]
                        t = 1.0/(abs(theta)+np.sqrt(1.0+theta**2))
                        if(theta < 0): t = -t
                    c = 1.0/np.sqrt(1+t**2)
                    s = t*c
                    tau = s/(1.0+c)
                    h = t*a[ip,iq]
                    z[ip] = z[ip] - h
                    z[iq] = z[iq] + h
                    d[ip] = d[ip] - h
                    d[iq] = d[iq] + h
                    a[ip,iq] = 0.0
                    for j in range(ip):
                        g = a[j,ip]
                        h = a[j,iq]
                        a[j,ip] = g - s*(h+g*tau)
                        a[j,iq] = h + s*(g-h*tau)
                    for j in range(ip+1,iq):
                        g = a[ip,j]
                        h = a[j,iq]
                        a[ip,j] = g - s*(h+g*tau)
                        a[j,iq] = h + s*(g-h*tau)
                    for j in range(iq+1,n):
                        g = a[ip,j]
                        h = a[iq,j]
                        a[ip,j] = g - s*(h+g*tau)
                        a[iq,j] = h + s*(g-h*tau)
                    for j in range(n):
                        g = v[j,ip]
                        h = v[j,iq]
                        v[j,ip] = g - s*(h+g*tau)
                        v[j,iq] = h + s*(g-h*tau)
                    nbinot = nbinot + 1
        for ip in range(n):
            b[ip] = b[ip] + z[ip]
            d[ip] = b[ip]
            z[ip] = 0.0
    raise RuntimeError(f"Too many iterations ({itermax}) in routine jacobi")

def jacobi_gpt1(a, nmax=500, itermax=50, copy=True):
    a = np.array(a, dtype=np.float64, copy=copy, order="C")
    n = a.shape[0]

    if a.ndim != 2 or n != a.shape[1]:
        raise ValueError("`a` must be a square 2D array")
    if n > nmax:
        raise BufferError(f"n(={n}) is too large (nmax={nmax})")

    b = np.diag(a).copy()
    d = b.copy()
    z = np.zeros(n, dtype=np.float64)
    v = np.eye(n, dtype=np.float64)

    nbinot = 0

    for i in range(itermax):
        # Same meaning as sum(abs(triu(a,1))) but without materializing triu(a,1)
        sm = np.abs(a[np.triu_indices(n, k=1)]).sum()
        if sm == 0.0:
            return d, v

        tresh = 0.2 * sm / (n * n) if i < 4 else 0.0

        for ip in range(n - 1):
            for iq in range(ip + 1, n):
                a_ipiq = a[ip, iq]
                g = 100.0 * abs(a_ipiq)

                if (i > 4) and ((abs(d[ip]) + g) == abs(d[ip])) and ((abs(d[iq]) + g) == abs(d[iq])):
                    a[ip, iq] = 0.0

                elif abs(a_ipiq) > tresh:
                    h = d[iq] - d[ip]

                    if (abs(h) + g) == abs(h):
                        t = a_ipiq / h
                    else:
                        theta = 0.5 * h / a_ipiq
                        t = 1.0 / (abs(theta) + np.sqrt(1.0 + theta * theta))
                        if theta < 0.0:
                            t = -t

                    c = 1.0 / np.sqrt(1.0 + t * t)
                    s = t * c
                    tau = s / (1.0 + c)
                    h = t * a_ipiq

                    z[ip] -= h
                    z[iq] += h
                    d[ip] -= h
                    d[iq] += h
                    a[ip, iq] = 0.0

                    # for j in range(ip):
                    if ip > 0:
                        gvec = a[:ip, ip].copy()
                        hvec = a[:ip, iq].copy()
                        a[:ip, ip] = gvec - s * (hvec + gvec * tau)
                        a[:ip, iq] = hvec + s * (gvec - hvec * tau)

                    # for j in range(ip+1, iq):
                    if iq - ip > 1:
                        gvec = a[ip, ip + 1:iq].copy()
                        hvec = a[ip + 1:iq, iq].copy()
                        a[ip, ip + 1:iq] = gvec - s * (hvec + gvec * tau)
                        a[ip + 1:iq, iq] = hvec + s * (gvec - hvec * tau)

                    # for j in range(iq+1, n):
                    if iq + 1 < n:
                        gvec = a[ip, iq + 1:n].copy()
                        hvec = a[iq, iq + 1:n].copy()
                        a[ip, iq + 1:n] = gvec - s * (hvec + gvec * tau)
                        a[iq, iq + 1:n] = hvec + s * (gvec - hvec * tau)

                    # for j in range(n):
                    gvec = v[:, ip].copy()
                    hvec = v[:, iq].copy()
                    v[:, ip] = gvec - s * (hvec + gvec * tau)
                    v[:, iq] = hvec + s * (gvec - hvec * tau)

                    nbinot += 1

        b += z
        d = b.copy()
        z.fill(0.0)

    raise RuntimeError(f"Too many iterations ({itermax}) in routine jacobi")


def jacobi_gemini1(a_in, itermax=50):
    a = a_in.copy().astype(float)
    n = a.shape[0]
    v = np.eye(n)
    
    for i in range(itermax):
        # 상삼각 성분 중 절댓값이 가장 큰 원소의 위치 탐색 (Pivot 선택)
        # 3x3이므로 직접 인덱싱이 효율적일 수 있으나 확장성을 위해 triu 유지
        upper_tri_indices = np.triu_indices(n, k=1)
        abs_upper = np.abs(a[upper_tri_indices])
        sm = np.sum(abs_upper)
        
        # 수렴 조건: 비대각 성분의 합이 충분히 작을 때
        if sm < 1e-15:
            return np.diag(a), v

        # 피벗 선택 (가장 큰 값을 가진 index 추출)
        idx = np.argmax(abs_upper)
        p, q = upper_tri_indices[0][idx], upper_tri_indices[1][idx]
        
        apq = a[p, q]
        if abs(apq) < 1e-15:
            continue

        # 회전각 계산
        h = a[q, q] - a[p, p]
        if abs(h) + abs(apq) * 100 == abs(h): # Underflow 방지
            t = apq / h
        else:
            theta = 0.5 * h / apq
            t = 1.0 / (np.abs(theta) + np.sqrt(1.0 + theta**2))
            if theta < 0: t = -t
            
        c = 1.0 / np.sqrt(1 + t**2)
        s = t * c
        tau = s / (1.0 + c)
        
        # 행렬 원소 업데이트 (벡터화)
        h = t * apq
        a[p, p] -= h
        a[q, q] += h
        a[p, q] = 0.0
        
        # p행/열과 q행/열의 다른 원소들 업데이트
        for j in range(n):
            if j != p and j != q:
                # 대칭 행렬 성질을 이용하되, 인덱스 관리를 위해 임시 저장
                # a[j, p]와 a[p, j]는 동일함 (대칭 가정)
                g = a[min(j, p), max(j, p)]
                h = a[min(j, q), max(j, q)]
                
                new_g = g - s * (h + g * tau)
                new_h = h + s * (g - h * tau)
                
                a[min(j, p), max(j, p)] = new_g
                a[min(j, q), max(j, q)] = new_h

        # 고유벡터 행렬 업데이트
        v_p = v[:, p].copy()
        v_q = v[:, q].copy()
        v[:, p] = v_p - s * (v_q + v_p * tau)
        v[:, q] = v_q + s * (v_p - v_q * tau)

    raise RuntimeError(f"최대 반복 횟수({itermax}) 내에 수렴하지 않았습니다.")

#***********************************************************************
# subroutine `indexx` can be replaced by `np.argsort`

#***********************************************************************
from scipy.special import elliprf
rf = elliprf
# def rf(x, y, z):
#     errtol = 0.08
#     tiny = 1.5e-38
#     big = 3.e37
#     third = 1. / 3.
#     c1 = 1. / 24.
#     c2 = 0.1
#     c3 = 3. / 44.
#     c4 = 1. / 14.

#     if min(x, y, z) < 0.0 or min(x + y, x + z, y + z) < tiny or max(x, y, z) > big:
#         raise ValueError("Invalid arguments in rf")

#     xt = x
#     yt = y
#     zt = z

#     while True:
#         sqrtx = np.sqrt(xt)
#         sqrty = np.sqrt(yt)
#         sqrtz = np.sqrt(zt)
#         alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz
#         xt = 0.25 * (xt + alamb)
#         yt = 0.25 * (yt + alamb)
#         zt = 0.25 * (zt + alamb)
#         ave = third * (xt + yt + zt)
#         delx = (ave - xt) / ave
#         dely = (ave - yt) / ave
#         delz = (ave - zt) / ave
#         if max(abs(delx), abs(dely), abs(delz)) > errtol:
#             continue
#         e2 = delx * dely - delz ** 2
#         e3 = delx * dely * delz
#         return (1. + (c1 * e2 - c2 - c3 * e3) * e2 + c4 * e3) / np.sqrt(ave)

#***********************************************************************
def cubic(a1, a2, a3, a4):
    """
    Function for evaluating a cubic equation. In the case where
    there is one real root, this routine returns it. If there are three
    roots, it returns the smallest positive one.

    The cubic equation is a1*x**3 + a2*x**2 + a3*x + a4 = 0.

    Ref: Tallarida, R. J., "Pocket Book of Integrals and Mathematical
    Formulas," 2nd ed., Boca Raton: CRC Press 1992, pp. 8-9.
    """

    if a1 == 0.0:
        raise ValueError('Quadratic/linear equation passed to cubic function')

    p = a3/a1 - (a2/a1)**2/3.0
    q = a4/a1 - a2*a3/a1**2/3.0 + 2.0*(a2/a1)**3/27.0
    d = p**3/27.0 + q**2/4.0

    if d > 0.0:
        a = -0.5*q + np.sqrt(d)
        b = -0.5*q - np.sqrt(d)
        a = a**(1/3) if a > 0 else -(-a)**(1/3)
        b = b**(1/3) if b > 0 else -(-b)**(1/3)
        y = a + b
    else:
        ar = -0.5 * q
        ai = np.sqrt(-d)
        r = (ar**2 + ai**2)**(1/6)
        theta = np.arctan2(ai, ar)
        y1 = 2.0 * r * np.cos(theta/3.0) - a2/a1/3.0
        y = y1
        y2 = 2.0 * r * np.cos(theta/3.0 + 2.0*np.pi/3.0) - a2/a1/3.0
        if y < 0 or (y2 > 0 and y2 < y): y = y2
        y3 = 2.0 * r * np.cos(theta/3.0 - 2.0*np.pi/3.0) - a2/a1/3.0
        if y < 0 or (y3 > 0 and y3 < y): y = y3

    return y

#***********************************************************************
# subroutine `int_indexx` can be replaced by `np.argsort`

#***********************************************************************
# function `ran2` can be replaced by `np.random.rand``


#=======================================================================
def spline(x):
#=======================================================================
#   Moved from `compute_neiKDtree_mod`
    if (x<=1.):
        ans=1.-1.5*x**2+0.75*x**3
    elif (x<=2.):
        ans=0.25*(2.-x)**3
    else:
        ans=0.
    return ans

#=======================================================================
def splines(xarr):
#=======================================================================
    ans = np.zeros_like(xarr)
    mask1 = (xarr <= 1.0)
    x1 = xarr[mask1]
    x12 = x1*x1
    mask2 = (~mask1) & (xarr <= 2.0)
    y = 2 - xarr[mask2]
    ans[mask1] = 1.0 - 1.5*x12 + 0.75*x12*x1
    ans[mask2] = 0.25 * y*y*y
    return ans

@njit(cache=True, fastmath=True, parallel=True)
def splines_numba(xarr):
    n, m = xarr.shape
    ans = np.empty_like(xarr)

    for ipar0 in prange(n):
        for i in range(m):
            x = xarr[ipar0, i]
            if x <= 1.0:
                x2 = x * x
                ans[ipar0, i] = 1.0 - 1.5 * x2 + 0.75 * x2 * x
            elif x <= 2.0:
                y = 2.0 - x
                y2 = y * y
                ans[ipar0, i] = 0.25 * y2 * y
            else:
                ans[ipar0, i] = 0.0

    return ans

#=======================================================================
def icellid(xtest:np.ndarray):
#=======================================================================
#  Compute cell id corresponding to the signs of coordinates of xtest
#  as follows :
#  (-,-,-) : 0
#  (+,-,-) : 1
#  (-,+,-) : 2
#  (+,+,-) : 3
#  (-,-,+) : 4
#  (+,-,+) : 5
#  (-,+,+) : 6
#  (+,+,+) : 7
#  For self-consistency, the array pos_ref_0 should be defined exactly 
#  with the same conventions
#=======================================================================
    # Moved from `compute_neiKDtree_mod`
    assert xtest.shape == (3,)
    icellid3d = (xtest>=0).astype(np.int32)
    return icellid3d[0]+2*icellid3d[1]+4*icellid3d[2]



# arr_for_icellids = 2 ** np.arange(3,dtype=np.int8)
arr_for_icellids = np.array([1,2,4], dtype=np.int8)
#=======================================================================
def icellids(xtests:np.ndarray, pos_this_node=None, mode=1):
    '''
    (N,3) -> (N,)
    '''
#=======================================================================
#  Compute cell id corresponding to the signs of coordinates of xtest
#  as follows :
#  (-,-,-) : 0
#  (+,-,-) : 1
#  (-,+,-) : 2
#  (+,+,-) : 3
#  (-,-,+) : 4
#  (+,-,+) : 5
#  (-,+,+) : 6
#  (+,+,+) : 7
#  For self-consistency, the array pos_ref_0 should be defined exactly 
#  with the same conventions
#=======================================================================
    global arr_for_icellids 
    if mode==0:
        # Moved from `compute_neiKDtree_mod`
        is_positive = xtests >= 0
        icellid = np.sum(is_positive * arr_for_icellids, axis=1, dtype=np.uint64)
        return icellid
    elif mode==1: # <--- this is the fastest
        if pos_this_node is None:
            b = (xtests >= 0).astype(np.uint64, copy=False)
        else:
            b = (xtests >= pos_this_node).astype(np.uint64, copy=False)
        return (b @ arr_for_icellids).astype(np.uint64, copy=False)

    elif mode==2:
        if pos_this_node is None:
            return (
                (xtests[:, 0] >= 0).astype(np.uint64)
                | ((xtests[:, 1] >= 0).astype(np.uint64) << 1)
                | ((xtests[:, 2] >= 0).astype(np.uint64) << 2)
            )
        else:
            return (
                (xtests[:, 0] >= pos_this_node[0]).astype(np.uint64)
                | ((xtests[:, 1] >= pos_this_node[1]).astype(np.uint64) << 1)
                | ((xtests[:, 2] >= pos_this_node[2]).astype(np.uint64) << 2)
            )


