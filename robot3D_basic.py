#!/usr/bin/env python
# coding: utf-8


from vedo import *
import numpy as np

def RotationMatrix(theta, axis_name):
    """ calculate single rotation of $theta$ matrix around x,y or z
        code from: https://programming-surgeon.com/en/euler-angle-python-en/
    input
        theta = rotation angle(degrees)
        axis_name = 'x', 'y' or 'z'
    output
        3x3 rotation matrix
    """

    c = np.cos(theta * np.pi / 180)
    s = np.sin(theta * np.pi / 180)
    
    if axis_name =='x':
        rotation_matrix = np.array([[1, 0,  0],
                                    [0, c, -s],
                                    [0, s,  c]])
    if axis_name =='y':
        rotation_matrix = np.array([[ c,  0, s],
                                    [ 0,  1, 0],
                                    [-s,  0, c]])
    elif axis_name =='z':
        rotation_matrix = np.array([[c, -s, 0],
                                    [s,  c, 0],
                                    [0,  0, 1]])
    return rotation_matrix



def createCoordinateFrameMesh():
    """Returns the mesh representing a coordinate frame
    Args:
      No input args
    Returns:
      F: vedo.mesh object (arrows for axis)
      
    """         
    _shaft_radius = 0.05
    _head_radius = 0.10
    _alpha = 1
    
    
    # x-axis as an arrow  
    x_axisArrow = Arrow(start_pt=(0, 0, 0),
                        end_pt=(1, 0, 0),
                        s=None,
                        shaft_radius=_shaft_radius,
                        head_radius=_head_radius,
                        head_length=None,
                        res=12,
                        c='red',
                        alpha=_alpha)

    # y-axis as an arrow  
    y_axisArrow = Arrow(start_pt=(0, 0, 0),
                        end_pt=(0, 1, 0),
                        s=None,
                        shaft_radius=_shaft_radius,
                        head_radius=_head_radius,
                        head_length=None,
                        res=12,
                        c='green',
                        alpha=_alpha)

    # z-axis as an arrow  
    z_axisArrow = Arrow(start_pt=(0, 0, 0),
                        end_pt=(0, 0, 1),
                        s=None,
                        shaft_radius=_shaft_radius,
                        head_radius=_head_radius,
                        head_length=None,
                        res=12,
                        c='blue',
                        alpha=_alpha)
    
    originDot = Sphere(pos=[0,0,0], 
                       c="black", 
                       r=0.10)


    # Combine the axes together to form a frame as a single mesh object 
    F = x_axisArrow + y_axisArrow + z_axisArrow + originDot
        
    return F


def getLocalFrameMatrix(R_ij, t_ij): 
    """Returns the matrix representing the local frame
    Args:
      R_ij: rotation of Frame j w.r.t. Frame i
      t_ij: translation of Frame j w.r.t. Frame i
    Returns:
      T_ij: Matrix of Frame j w.r.t. Frame i. 
      
    """             
    # Rigid-body transformation [ R t ]
    T_ij = np.block([[R_ij,                t_ij],
                     [np.zeros((1, 3)),       1]])
    
    return T_ij
   

def forward_kinematics(phi, L1, L2, L3, L4):

 # Matrix of Frame 1 (written w.r.t. Frame 0, which is the previous frame) 
    R_01 = RotationMatrix(phi[0], axis_name = 'z')   # Rotation matrix
    p1   = np.array([[3],[2], [0.0]])              # Frame's origin (w.r.t. previous frame)
    t_01 = p1                                      # Translation vector

    T_01 = getLocalFrameMatrix(R_01, t_01)         # Matrix of Frame 1 w.r.t. Frame 0 (i.e., the world frame)


    # Matrix of Frame 2 (written w.r.t. Frame 1, which is the previous frame)   
    R_12 = RotationMatrix(phi[1], axis_name = 'z')   # Rotation matrix
    p2   = np.array([[L1+0.8], [0.0], [0.0]])           # Frame's origin (w.r.t. previous frame)
    t_12 = p2                                      # Translation vector

    # Matrix of Frame 2 w.r.t. Frame 1 
    T_12 = getLocalFrameMatrix(R_12, t_12)

    # Matrix of Frame 2 w.r.t. Frame 0 (i.e., the world frame)
    T_02 = T_01 @ T_12



    # Matrix of Frame 3 (written w.r.t. Frame 2, which is the previous frame)   
    R_23 = RotationMatrix(phi[2], axis_name = 'z')   # Rotation matrix
    p3   = np.array([[L2+0.8],[0.0], [0.0]])           # Frame's origin (w.r.t. previous frame)
    t_23 = p3                                      # Translation vector

    # Matrix of Frame 3 w.r.t. Frame 2 
    T_23 = getLocalFrameMatrix(R_23, t_23)

    # Matrix of Frame 3 w.r.t. Frame 0 (i.e., the world frame)
    T_03 = T_02 @ T_23


    R_34 = RotationMatrix(phi[3], axis_name = 'z')   # Rotation matrix
    p4   = np.array([[L3+0.4],[0.0], [0.0]])           # Frame's origin (w.r.t. previous frame)
    t_34 = p4                                      # Translation vector

    # Matrix of Frame 3 w.r.t. Frame 2 
    T_34 = getLocalFrameMatrix(R_34, t_34)

    # Matrix of Frame 3 w.r.t. Frame 0 (i.e., the world frame)
    T_04 = T_03 @ T_34



    return T_01, T_02, T_03, T_04, T_04[:-1,-1]


class Arm:
    def __init__(self, ax, f1, f2, f3, f4):
        self.ax = ax
        self.f1 = f1
        self.f2 = f2
        self.f3 = f3
        self.f4 = f4

    def unpack(self):
        return [self.ax, self.f1, self.f2, self.f3, self.f4]

def vid_runner(i):
    phis = np.vstack([np.linspace(0, 30, 50), np.linspace(0, -50, 50), np.linspace(0, -30, 50), np.zeros(50)]).T
    # Set the limits of the graph x, y, and z ranges 
    axes = Axes(xrange=(0,20), yrange=(-2,10), zrange=(0,6))

    # Lengths of arm parts 
    L1 = 5   # Length of link 1
    L2 = 8   # Length of link 2
    L3 = 3
    L4 = 0

    # Joint angles 

    T_01, T_02, T_03, T_04, e = forward_kinematics(phis[i], L1, L2, L3, L4)
    # Create the coordinate frame mesh and transform
    Frame1Arrows = createCoordinateFrameMesh()
    
    # Now, let's create a cylinder and add it to the local coordinate frame
    r1 = 0.4
    link1_mesh = Cylinder(r=0.4, 
                          height=L1, 
                          pos = (r1+L1/2,0,0),
                          c="blue", 
                          alpha=.8, 
                          axis=(1,0,0)
                          )
    
    # Also create a sphere to show as an example of a joint
    sphere1 = Sphere(r=r1).pos(0,0,0).color("gray").alpha(.8)

    # Combine all parts into a single object 
    Frame1Arrows.apply_transform(T_01)  
    link1_mesh.apply_transform(T_01)  
    sphere1.apply_transform(T_01)  
    Frame1 = Group((Frame1Arrows, link1_mesh, sphere1))

    # Transform the part to position it at its correct location and orientation 
    
    # Create the coordinate frame mesh and transform
    Frame2Arrows = createCoordinateFrameMesh()
    sphere2 = Sphere(r=r1).pos(0,0,0).color("gray").alpha(.8)
    # Now, let's create a cylinder and add it to the local coordinate frame
    link2_mesh = Cylinder(r=0.4, 
                          height=L2, 
                          pos = (r1+L2/2,0,0),
                          c="blue", 
                          alpha=.8, 
                          axis=(1,0,0)
                          )

    # Combine all parts into a single object 
    Frame2Arrows.apply_transform(T_02)
    link2_mesh.apply_transform(T_02)
    sphere2.apply_transform(T_02)

    Frame2 = Group((Frame2Arrows, link2_mesh, sphere2))

    # Transform the part to position it at its correct location and orientation 
    # Create the coordinate frame mesh and transform. This point is the end-effector. So, I am 
    # just creating the coordinate frame. 
    Frame3Arrows = createCoordinateFrameMesh()
    sphere3 = Sphere(r=r1).pos(0,0,0).color("gray").alpha(.8)

    # Now, let's create a cylinder and add it to the local coordinate frame
    link3_mesh = Cylinder(r=0.4, 
                          height=L3, 
                          pos = (r1+L3/2,0,0),
                          c="blue", 
                          alpha=.8, 
                          axis=(1,0,0)
                          )
    Frame3Arrows.apply_transform(T_03)
    link3_mesh.apply_transform(T_03)
    sphere3.apply_transform(T_03)
    # Combine all parts into a single object 
    Frame3 = Group((Frame3Arrows, link3_mesh, sphere3))

    # Transform the part to position it at its correct location and orientation 
    # Create the coordinate frame mesh and transform. This point is the end-effector. So, I am 
    # just creating the coordinate frame. 
    Frame4Arrows = createCoordinateFrameMesh()

    # Transform the part to position it at its correct location and orientation 
    Frame4Arrows.apply_transform(T_04) 
    Frame4 = Frame4Arrows
    
    plt = Plotter()
    plt += [axes, Frame1, Frame2, Frame3, Frame4]
    plt.render()
    plt.show(screenshot="".join(("pics/",str(i),".png"))).close()
    
    # Show everything 
    # show([Frame1, Frame2, Frame3, Frame4], axes, viewup="z").close()
    # show([Frame1Arrows], axes, viewup="z").close()




def main():

    # Set the limits of the graph x, y, and z ranges 
    axes = Axes(xrange=(0,20), yrange=(-2,10), zrange=(0,6))

    # Lengths of arm parts 
    L1 = 5   # Length of link 1
    L2 = 8   # Length of link 2
    L3 = 3
    L4 = 0

    # Joint angles 
    phi = np.array([30,-50,-30,0])

    T_01, T_02, T_03, T_04, e = forward_kinematics(phi, L1, L2, L3, L4)
    # Create the coordinate frame mesh and transform
    Frame1Arrows = createCoordinateFrameMesh()
    
    # Now, let's create a cylinder and add it to the local coordinate frame
    r1 = 0.4
    link1_mesh = Cylinder(r=0.4, 
                          height=L1, 
                          pos = (r1+L1/2,0,0),
                          c="yellow", 
                          alpha=.8, 
                          axis=(1,0,0)
                          )
    
    # Also create a sphere to show as an example of a joint
    sphere1 = Sphere(r=r1).pos(0,0,0).color("gray").alpha(.8)

    # Combine all parts into a single object 
    Frame1Arrows.apply_transform(T_01)  
    link1_mesh.apply_transform(T_01)  
    sphere1.apply_transform(T_01)  
    Frame1 = Group((Frame1Arrows, link1_mesh, sphere1))

    # Transform the part to position it at its correct location and orientation 
    
    # Create the coordinate frame mesh and transform
    Frame2Arrows = createCoordinateFrameMesh()
    sphere2 = Sphere(r=r1).pos(0,0,0).color("gray").alpha(.8)
    # Now, let's create a cylinder and add it to the local coordinate frame
    link2_mesh = Cylinder(r=0.4, 
                          height=L2, 
                          pos = (r1+L2/2,0,0),
                          c="red", 
                          alpha=.8, 
                          axis=(1,0,0)
                          )

    # Combine all parts into a single object 
    Frame2Arrows.apply_transform(T_02)
    link2_mesh.apply_transform(T_02)
    sphere2.apply_transform(T_02)

    Frame2 = Group((Frame2Arrows, link2_mesh, sphere2))

    # Transform the part to position it at its correct location and orientation 
    # Create the coordinate frame mesh and transform. This point is the end-effector. So, I am 
    # just creating the coordinate frame. 
    Frame3Arrows = createCoordinateFrameMesh()
    sphere3 = Sphere(r=r1).pos(0,0,0).color("gray").alpha(.8)

    # Now, let's create a cylinder and add it to the local coordinate frame
    link3_mesh = Cylinder(r=0.4, 
                          height=L3, 
                          pos = (r1+L3/2,0,0),
                          c="red", 
                          alpha=.8, 
                          axis=(1,0,0)
                          )
    Frame3Arrows.apply_transform(T_03)
    link3_mesh.apply_transform(T_03)
    sphere3.apply_transform(T_03)
    # Combine all parts into a single object 
    Frame3 = Group((Frame3Arrows, link3_mesh, sphere3))

    # Transform the part to position it at its correct location and orientation 
    # Create the coordinate frame mesh and transform. This point is the end-effector. So, I am 
    # just creating the coordinate frame. 
    Frame4Arrows = createCoordinateFrameMesh()

    # Transform the part to position it at its correct location and orientation 
    Frame4Arrows.apply_transform(T_04) 
    Frame4 = Frame4Arrows
    plt = Plotter()

    plt += [axes, Frame1, Frame2, Frame3, Frame4]
    plt.render()
    plt.show().close()
    # Show everything 
    # show([Frame1, Frame2, Frame3, Frame4], axes, viewup="z").close()
    # show([Frame1Arrows], axes, viewup="z").close()
    


if __name__ == '__main__':
    main()



