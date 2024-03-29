import math

def deltaR2( e1, p1, e2, p2):

    de = e1 - e2
    dp = deltaPhi(p1, p2)
    return de*de + dp*dp
            
def deltaR( *args ):

    return math.sqrt( deltaR2(*args) )
                
def deltaPhi( p1, p2):

    res = p1 - p2
    while res > math.pi:
        res -= 2*math.pi
    while res < -math.pi:
        res += 2*math.pi
    return res
