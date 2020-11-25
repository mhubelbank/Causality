import Behavior

print('Pa_B = ', .95, ', Pa_nB = ', .05)

def PB(K):
    Pa = [1] * K + [0] * (5-K)
    Pa_B = [.95] * 5
    Pa_nB = [.05] * 5
    PB = .005
    print('Pa = ', Pa, ', Pa_B = ', Pa_B, ', Pa_nB = ', Pa_nB, ', PB = ', PB)
    result = Behavior.calcProb(Pa, Pa_B, Pa_nB, PB)
    return result


for k in range(1, 6):
    pb = PB(k)
    print('K = ', k, 'P(B) = ', pb)
