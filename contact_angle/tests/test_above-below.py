import contact_angle as cnt
import mdtraj as md
import numpy as np

from contact_angle.utils.general import get_fn


def test_flipped():
    """Same trajectory, but rotated around x-axis - should give same contact angle
    """
    traj = md.load(get_fn('chol-tail-original.dcd'),
            top=get_fn('chol-wetting.hoomdxml'))
    ca = cnt.calc_contact_angle(traj.xyz[:, 14400:]*6, guess_R=4.0, guess_z0=1.8, 
            guess_rho_n=1.0, left_tol=0.1, z_range=(-0.1, 9), surface_normal='z', 
            n_bins=100, fit_range=(2, 4.0), droplet_location='above')
    traj = md.load(get_fn('chol-tail-rotated.dcd'),
            top=get_fn('chol-wetting.hoomdxml'))
    ca2 = cnt.calc_contact_angle(traj.xyz[:, 14400:]*6, guess_R=4.0, guess_z0=1.8, 
            guess_rho_n=1.0, left_tol=0.1, z_range=(-7.0, 0.3), surface_normal='z', 
            n_bins=100, fit_range=(-6.0, -1.0), droplet_location='below')
    assert(np.absolute(ca['theta'] - ca2['theta']) < 5.0)
