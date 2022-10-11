"""
Functions describing flutter
"""

def flutter(flutter_velocity, a_speed_sound, ambient_pressure, ground_pressure, chord_root, Aspect_ratio, G_shear, chord_ratio):
    """
    Returns the thickness of a planar square sectioned fin such that it has a flutter velocity equal to a pre-defined value.
    Based on NACA 4197 fin flutter. Fin flutter is worst at Max Q.

    Parameters
    ----------
    flutter_velocity : float
      velocity up to which the fins are safe from flutter
    a_speed_sound : float
      velocity of the speed of sound
    ambient_pressure : float
        external pressure at maximum Q, static
    ground_pressure : float
        ground/launch site pressure, static
    chord_root : float
        root chord of the fin
    G_shear : float
        shear modulus of the fin material in psi
    chord_ratio: float
        ratio of Chord_tip/Chord_root
      
    Returns
    -------
    Cp_x_distance : float
      distance between the top of the fin and centre of pressure for supersonic regime
    """
    result = (((flutter_velocity/a_speed_sound)**2)*39.3*(ambient_pressure/ground_pressure)*((chord_ratio+1)/2)*(Aspect_ratio**3))/(G_shear*(Aspect_ratio+2))
    t = chord_root * ((result)**(1./3))
    return t