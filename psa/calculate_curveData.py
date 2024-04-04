import numpy as np

def number_of_pulses(PS_mat, inpoFact):
    number_of_pulses = PS_mat.shape[0]/inpoFact
    return number_of_pulses 

def angle_between_vectors(vector1, vector2):
    # Calculate dot product of the two vectors
    dot_product = np.dot(vector1, vector2)
    
    # Calculate magnitudes of the vectors
    magnitude1 = np.linalg.norm(vector1)
    magnitude2 = np.linalg.norm(vector2)
    
    # Calculate the angle between the vectors using the dot product formula
    angle_rad = np.arccos(dot_product / (magnitude1 * magnitude2))
    
    # Convert angle from radians to degrees
    angle_deg = np.degrees(angle_rad)
    
    return angle_deg
    
def angle_with_x_axis(vector):
    # Project the vector onto the xy-plane
    projected_vector = np.array([vector[0], vector[1], 0])
    
    # Calculate the angle between the projected vector and the x-axis
    angle_rad = np.arctan2(projected_vector[1], projected_vector[0])
    
    # Convert angle from radians to degrees
    angle_deg = np.degrees(angle_rad)
    
    return angle_deg
    
	
def calculate_closed_curve_area_app(CM2D, close_curve=False):
    """
    calculate_closed_curve_area_app calculates the area of a 2-dimensional curve
    :param CM2D: CM of a 2-dimensional curve (projection of original CM on a plane)
    :param close_curve: if Area should be calculated from unchanged curve or from a closed version of the curve (last point = first point)
    :return: Area of the curve
    """
    if close_curve:
        CM2D = np.vstack((CM2D, CM2D[0]))  # closes curve if close_curve is True

    A = np.trapz(CM2D[:, 1], x=CM2D[:, 0])  # calculate the area using trapezoidal integration
    return A
    
def curvature_torsion_3d(curve):
    # Calculate differences between successive points
    dp = np.diff(curve, axis=0)
    dp_length = np.linalg.norm(dp, axis=1)
    
    # Calculate derivatives of differences
    dp1 = np.gradient(dp, axis=0)
    dp1_length = np.linalg.norm(dp1, axis=1)
    
    # Calculate curvature (using norm of second derivative)
    curvature = np.divide(dp1_length, dp_length**2)
    
    # Calculate torsion
    dp2 = np.gradient(dp1, axis=0)
    torsion = np.sum(dp * np.cross(dp1, dp2), axis=1) / (dp_length**2 * dp1_length**2)
    
    # Calculate the arc length of the curve
    arc_length = np.cumsum(dp_length)

    # Calculate integrated curvature and torsion
    integrated_curvature = np.trapz(curvature, x=arc_length)
    integrated_torsion = np.trapz(torsion, x=arc_length)
    integrated_absolut_torsion = np.trapz(abs(torsion), x=arc_length)
    
    # Calculate average curvature and torsion
    avg_curvature = np.mean(curvature)
    avg_torsion = np.mean(torsion)
    
    return arc_length, curvature, torsion, integrated_curvature, integrated_torsion, integrated_absolut_torsion, avg_curvature, avg_torsion

