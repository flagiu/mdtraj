import numpy as np

def read_sqt(filename: str, dt=0.002):
    """
    Read wavevectors q, time_intervals t, and S(q,t) and its error from file.
    - dt: timestep
    """
    f=open(filename,"r")
    lines=f.readlines()
    # Parse wavevectors from 1st block:
    i = 0
    block=0
    q = []
    while block==0:
        line = lines[i]
        words = line.strip('\n').split(' ')
        if line[0]=='#':
            words=None
        else:
            assert len(words)==1
            if words[0]=='':
                block=1 # separation empty line
            else:
                q.append(float(words[0]))
        i+=1
    q = np.array(q)
    Nq=len(q)

    # Parse delta timesteps from 2nd block:
    t = []
    while block==1:
        line = lines[i]
        words = line.strip('\n').split(' ')
        if line[0]=='#':
            words=None
        else:
            assert len(words)==1
            if words[0]=='':
                block=2 # separation empty line
            else:
                t.append(float(words[0]))
        i+=1
    t = dt * np.array(t)
    Nt=len(t)
    f.close()

    # Load 3rd block
    sqt = np.genfromtxt(filename, skip_header=1+Nq+1+Nt+1, skip_footer=Nq)
    # Load 4th block
    sqt_ = np.genfromtxt(filename, skip_header=1+Nq+1+Nt+1+Nq+1, skip_footer=0)
    if Nq==1 or Nt==1:
        assert len(sqt)==Nt or len(sqt)==Nq
        assert len(sqt_)==Nt or len(sqt_)==Nq
        sqt = sqt.reshape(Nq,Nt)
        sqt_ = sqt_.reshape(Nq,Nt)
    else:
        assert sqt.shape==(Nq,Nt)
        assert sqt_.shape==(Nq,Nt)

    return q,t,sqt,sqt_
