"""
Definition of the variational form of the motor problem
"""
from fea_utils import *
from piecewise_permeability import *


# START NEW PERMEABILITY
def RelativePermeability(subdomain, u, uhat):
    gradu = gradx(u,uhat)
    # if subdomain == 1: # Electrical/Silicon/Laminated Steel
    if subdomain == 1 or subdomain == 2: # Electrical/Silicon/Laminated Steel
        B = as_vector((gradu[1], -gradu[0]))
        norm_B = sqrt(dot(B, B) + DOLFIN_EPS)

        mu = conditional(
            lt(norm_B, 1.004),
            linearPortion(norm_B),
            conditional(
                lt(norm_B, 1.433),
                cubicPortion(norm_B),
                (expA * exp(expB*norm_B + expC) + 1)
            )
        )
        # mu = 4000. # constant case
    elif subdomain == 3:
        mu = 1.00 # insert value for titanium or shaft material
    elif subdomain >= 4 and subdomain <= 28: # AIR
        mu = 1.0
    elif subdomain >= 29 and subdomain <= 40: # NEODYMIUM
        mu = 1.05
    elif subdomain >= 41: # COPPER
        mu = 1.00

    return mu
# END NEW PERMEABILITY

def compute_i_abc(iq, angle=0.0):
    i_abc = as_vector([
        iq * np.sin(angle),
        iq * np.sin(angle + 2*np.pi/3),
        iq * np.sin(angle - 2*np.pi/3),
    ])
    return i_abc
    
def JS(v,uhat,iq,p,s,Hc,angle):
    """
    The variational form for the source term (current) of the
    Maxwell equation
    """
    Jm = 0.
    gradv = gradx(v,uhat)
    base_magnet_dir = 2 * np.pi / p / 2
    magnet_sweep    = 2 * np.pi / p
    for i in range(p):
        angle = base_magnet_dir + i * magnet_sweep
        Hx = Constant((-1)**i * Hc * np.cos(angle))
        Hy = Constant((-1)**i * Hc * np.sin(angle))

        H = as_vector([Hx, Hy])

        curl_v = as_vector([gradv[1],-gradv[0]])
        Jm += inner(H,curl_v)*dx(i + 4 + p*2 + 1)

    num_windings = 2*s
    num_phases = 3
    coil_per_phase = 6
    stator_winding_index_start  = 4 + 3 * p + 1
    stator_winding_index_end    = stator_winding_index_start + num_windings
    Jw = 0.
    N = 13
    i_abc = compute_i_abc(iq, angle)
    JA, JB, JC = i_abc[0] * N + DOLFIN_EPS, i_abc[1] * N + DOLFIN_EPS, i_abc[2] * N + DOLFIN_EPS

    # OLD METHOD
    # for i in range(int((num_windings) / (num_phases * coil_per_phase))):
    #     coil_start_ind = i * num_phases * coil_per_phase
    #     for j in range(3):
    #         if j == 0:
    #             J = JA
    #         elif j == 1:
    #             J = JB
    #         elif j == 2:
    #             J = JC
    #         phase_start_ind = coil_per_phase * j
    #         J_list = [
    #             J * (-1)**i * v * dx(stator_winding_index_start + phase_start_ind + coil_start_ind),
    #             J * (-1)**(i+1) * v * dx(stator_winding_index_start + 1 + phase_start_ind + coil_start_ind),
    #             J * (-1)**(i+1) * v * dx(stator_winding_index_start + 2 + phase_start_ind + coil_start_ind),
    #             J * (-1)**i * v * dx(stator_winding_index_start + 3 + phase_start_ind + coil_start_ind),
    #             J * (-1)**i * v * dx(stator_winding_index_start + 4 + phase_start_ind + coil_start_ind),
    #             J * (-1)**(i+1) * v * dx(stator_winding_index_start + 5 + phase_start_ind + coil_start_ind),
    #         ]
    #         Jw += sum(J_list)
            # order: + - - + + - (signs switch with each instance of the phases)

    # NEW METHOD
    coil_per_phase = 2
    num_windings = s
    for i in range(int((num_windings) / (num_phases * coil_per_phase))):
        coil_start_ind = i * num_phases * coil_per_phase
        for j in range(3):
            if i%3 == 0: # PHASE A
                J = JA
            if i%3 == 1: # PHASE C
                J = JB
            if i%3 == 2: # PHASE B
                J = JC
        
        J_list = [
            JB * (-1)**(i+1) * v * dx(stator_winding_index_start + coil_start_ind),
            JA * (-1)**(i) * v * dx(stator_winding_index_start + coil_start_ind + 1),
            JC * (-1)**(i+1) * v * dx(stator_winding_index_start + coil_start_ind + 2),
            JB * (-1)**(i) * v * dx(stator_winding_index_start + coil_start_ind + 3),
            JA * (-1)**(i+1) * v * dx(stator_winding_index_start + coil_start_ind + 4),
            JC * (-1)**(i) * v * dx(stator_winding_index_start + coil_start_ind + 5)
        ]
        Jw += sum(J_list)

    return Jm + Jw

def pdeRes(u,v,uhat,iq,dx,p,s,Hc,vacuum_perm,angle):
    """
    The variational form of the PDE residual for the magnetostatic problem
    """
    res = 0.
    gradu = gradx(u,uhat)
    gradv = gradx(v,uhat)
    num_components = 4 * 3 * p + 2 * s
    for i in range(num_components):
        res += 1./vacuum_perm*(1/RelativePermeability(i + 1, u, uhat))\
                *dot(gradu,gradv)*J(uhat)*dx(i + 1)
    res -= JS(v,uhat,iq,p,s,Hc,angle)
    return res

