package PythagoraDevice
  import SI = Modelica.SIunits;
  import Modelica.Math.Vectors.*;
  import Cv = Modelica.SIunits.Conversions;
  import C = Modelica.Constants;
  import Modelica.Mechanics.MultiBody.Frames;
  import Modelica.Mechanics.MultiBody.Types;
  import Modelica.Mechanics.MultiBody.*;
  import Modelica.Constants.*;

  package Visualizers
    extends Modelica.Icons.Package;

    model Sphere
      import Modelica.Mechanics.MultiBody.Types;
      import Modelica.Mechanics.MultiBody.Interfaces;
      import SI = Modelica.SIunits;
      import Modelica.Mechanics.MultiBody.Frames;
      Frames.Orientation R;
      parameter SI.Diameter sphereDiameter;
      parameter Types.Color sphereColor = {255, 255, 204};
      parameter String shapeType = "sphere";
      parameter Real r_shape[3];
      //  Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape sphere(shapeType = shapeType, color = sphereColor, length = sphereDiameter, width = sphereDiameter, height = sphereDiameter, lengthDirection = {1, 0, 0}, widthDirection = {0, 1, 0}, r_shape = r_shape, r = r_0, R = R);
//      Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape sphere(shapeType = "file://C:\Work\2020\Pita\DXF_Study\Ball\q78u2sv9fx1c-soccerball (1)\untitled.stl", color = sphereColor, length = sphereDiameter, width = sphereDiameter, height = sphereDiameter, lengthDirection = {1, 0, 0}, widthDirection = {0, 1, 0}, r_shape = {0, 0, 0} * sphereDiameter / 2, r = r_0, R = R);
            Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape sphere(shapeType = "sphere", color = sphereColor, length = sphereDiameter, width = sphereDiameter, height = sphereDiameter, lengthDirection = {1, 0, 0}, widthDirection = {0, 1, 0}, r_shape = {0, 0, 0} * sphereDiameter / 2, r = r_0, R = R);
      //    Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape sphere(shapeType="file://C:\Work\2020\Pita\DXF_Study\Ball\untitled.stl", color = sphereColor, length = sphereDiameter, width = sphereDiameter, height = sphereDiameter, lengthDirection = {1, 0, 0}, widthDirection = {0, 1, 0}, r_shape = -{1, 0, 0} * sphereDiameter / 2, r = r_0, R = R);
      parameter StateSelect stateSelect = StateSelect.avoid;
      SI.Position r_0[3](start = {0, 0, 0}, each stateSelect = stateSelect);
      PythagoraDevice.Interfaces.ShapeInput shapeInput annotation(
        Placement(visible = true, transformation(origin = {-106, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-110, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      r_0 = shapeInput.xyz;
      R = shapeInput.R;
    protected
      annotation(
        uses(Modelica(version = "3.2.3")),
        Icon(graphics = {Ellipse(fillColor = {0, 170, 255}, fillPattern = FillPattern.Sphere, extent = {{-100, 100}, {100, -100}}, endAngle = 360)}, coordinateSystem(initialScale = 0.1)));
    end Sphere;

    model Box
      import Modelica.Mechanics.MultiBody.Types;
      import Modelica.Mechanics.MultiBody.Interfaces;
      import SI = Modelica.SIunits;
      import Modelica.Mechanics.MultiBody.Frames;
      parameter Types.Color sphereColor = {255, 255, 255};
      parameter SI.Length Length = 1;
      parameter SI.Length width = 0.01;
      parameter SI.Length height = 0.01;
      parameter Real specularCoefficient = 1;
      parameter Real lengthDirection[3] = {1, 0, 0};
      parameter Real widthDirection[3] = {0, 1, 0};
      parameter Real r_shape[3] = {0, 0, 0};
      parameter Frames.Orientation R;
      Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape box(shapeType = "box", color = sphereColor, specularCoefficient = specularCoefficient, length = Length, width = width, height = height, lengthDirection = lengthDirection, widthDirection = widthDirection, r_shape = r_shape, r = r_0, R = R);
      Modelica.Blocks.Interfaces.RealInput r_0[3] annotation(
        Placement(visible = true, transformation(origin = {-114, 6}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-114, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    equation

      annotation(
        uses(Modelica(version = "3.2.3")),
        Icon(coordinateSystem(initialScale = 0.1), graphics = {Rectangle(fillColor = {255, 85, 127}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-100, 20}, {100, -20}})}));
    end Box;
  end Visualizers;

  package Examples
    extends Modelica.Icons.ExamplesPackage;

    model BallInBox "Test model"
      parameter Integer n = 45 "number of moving Ball";
      parameter Integer n_fix = 4 "number of fix ball";
      parameter Integer m = 3 "number of Wall";
      PythagoraDevice.Components.MultiBalls moveBalls(D = fill(0.1, n), m = fill(1, n), n = n, position = {{-0.4, 0.9, 0}, {-0.3, 0.9, 0}, {-0.2, 0.9, 0}, {-0.1, 0.9, 0}, {0, 0.9, 0}, {0.1, 0.9, 0}, {0.2, 0.9, 0}, {0.3, 0.9, 0}, {0.4, 0.9, 0}, {-0.4, 1.1, 0}, {-0.3, 1.1, 0}, {-0.2, 1.1, 0}, {-0.1, 1.1, 0}, {0, 1.1, 0}, {0.1, 1.1, 0}, {0.2, 1.1, 0}, {0.3, 1.1, 0}, {0.4, 1.1, 0}, {-0.4, 1.3, 0}, {-0.3, 1.3, 0}, {-0.2, 1.3, 0}, {-0.1, 1.3, 0}, {0, 1.3, 0}, {0.1, 1.3, 0}, {0.2, 1.3, 0}, {0.3, 1.3, 0}, {0.4, 1.3, 0}, {-0.4, 1.5, 0}, {-0.3, 1.5, 0}, {-0.2, 1.5, 0}, {-0.1, 1.5, 0}, {0, 1.5, 0}, {0.1, 1.5, 0}, {0.2, 1.5, 0}, {0.3, 1.5, 0}, {0.4, 1.5, 0}, {-0.4, 1.7, 0}, {-0.3, 1.7, 0}, {-0.2, 1.7, 0}, {-0.1, 1.7, 0}, {0, 1.7, 0}, {0.1, 1.7, 0}, {0.2, 1.7, 0}, {0.3, 1.7, 0}, {0.4, 1.7, 0}}, v_start = fill(0, n, 3)) annotation(
        Placement(visible = true, transformation(origin = {0, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.Contacts contacts[n * m] annotation(
        Placement(visible = true, transformation(origin = {26, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.FixedaWall fixedaWall[m](height = {1, 1, 1}, length = {1, 0.01, 0.01}, s0 = {{0, 0, 0}, {0.5, 1, 0}, {-0.5, 1, 0}}, width = {0.01, 2, 2}) annotation(
        Placement(visible = true, transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.FixedMultiBalls fixBalls(D = fill(0.07, n_fix), m = fill(1, n_fix), n = n_fix, position = {{-0.35, 0.6, 0}, {-0.15, 0.6, 0}, {0.05, 0.6, 0}, {0.25, 0.6, 0}}, v_start = fill(0, n_fix, 3)) annotation(
        Placement(visible = true, transformation(origin = {0, -18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.Contacts contactsFix[n * n_fix] annotation(
        Placement(visible = true, transformation(origin = {24, -18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
//Wall-movingBall
      for i in 1:n loop
        for j in 1:m loop
          connect(moveBalls.Center_a[i], contacts[(i - 1) * m + j].Center_b);
          connect(fixedaWall[j].Center_b, contacts[(i - 1) * m + j].Center_a);
        end for;
      end for;
//fixBall-movingBall
      for k in 1:n loop
        for p in 1:n_fix loop
          connect(moveBalls.Center_a[k], contactsFix[(k - 1) * n_fix + p].Center_b);
          connect(fixBalls.Center_a[p], contactsFix[(k - 1) * n_fix + p].Center_a);
        end for;
      end for;
      annotation(
        experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-06, Interval = 0.04),
        __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian -d=bltdump -d=bltdump -d=bltdump ",
        __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "rungekutta"));
    end BallInBox;
    
    model Contact_BallBall
      PythagoraDevice.Components.Ball ball1(g_f = {0, 0, 0}, v_g(fixed = true), v_g_ini = {5e-3, 0, 0}, w_g_ini = {0, 0, 0}) annotation(
        Placement(visible = true, transformation(origin = {-68, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.ContactsBallBall contactsAddRotation annotation(
        Placement(visible = true, transformation(origin = {-38, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.Ball ball2(g_f = {0, 0, 0}, position_ini = {0.015, 0, 0}, v_g(fixed = true), v_g_ini = {0, 0, 0}) annotation(
        Placement(visible = true, transformation(origin = {-4, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.ContactsBallBall contactsBallBall annotation(
        Placement(visible = true, transformation(origin = {-38, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.Ball ball4(g_f = {0, 0, 0}, position_ini = {0.015, -0.02, 0}, v_g(fixed = true), v_g_ini = {0, 0, 0}) annotation(
        Placement(visible = true, transformation(origin = {-4, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.Ball ball3(g_f = {0, 0, 0}, position_ini = {0, -0.02, 0}, v_g(fixed = true), v_g_ini = {5e-3, 0, 0}, w_g_ini = {0, 0, -5}) annotation(
        Placement(visible = true, transformation(origin = {-68, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.Ball ball5(g_f = {0, 0, 0}, position_ini = {0, -0.04, 0}, v_g(fixed = true), v_g_ini = {5e-3, 0, 0}, w_g_ini = {0, 0, -5}) annotation(
        Placement(visible = true, transformation(origin = {-68, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.Ball ball6(g_f = {0, 0, 0}, position_ini = {0.015, -0.04, 0}, v_g(fixed = true), v_g_ini = {0, 0, 0}, w_g_ini = {0, 0, -5}) annotation(
        Placement(visible = true, transformation(origin = {-4, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.ContactsBallBall contactsBallBall1 annotation(
        Placement(visible = true, transformation(origin = {-38, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.ContactsBallBall contactsBallBall2 annotation(
        Placement(visible = true, transformation(origin = {-38, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.Ball ball8(g_f = {0, 0, 0}, position_ini = {0.015, -0.06, 0}, v_g(fixed = true), v_g_ini = {0, 0, 0}, w_g_ini = {0, 0, 5}) annotation(
        Placement(visible = true, transformation(origin = {-4, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.Ball ball7(g_f = {0, 0, 0}, position_ini = {0, -0.06, 0}, v_g(fixed = true), v_g_ini = {5e-3, 0, 0}, w_g_ini = {0, 0, -5}) annotation(
        Placement(visible = true, transformation(origin = {-68, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(ball1.center, contactsAddRotation.center_a) annotation(
        Line(points = {{-68, 60}, {-48, 60}}, color = {0, 127, 0}));
      connect(contactsAddRotation.center_b, ball2.center) annotation(
        Line(points = {{-28, 60}, {-4, 60}}, color = {0, 79, 0}));
      connect(ball3.center, contactsBallBall.center_a) annotation(
        Line(points = {{-68, 20}, {-48, 20}}, color = {0, 127, 0}));
      connect(contactsBallBall.center_b, ball4.center) annotation(
        Line(points = {{-28, 20}, {-4, 20}}, color = {0, 79, 0}));
      connect(ball5.center, contactsBallBall1.center_a) annotation(
        Line(points = {{-68, -20}, {-48, -20}}, color = {0, 127, 0}));
      connect(contactsBallBall1.center_b, ball6.center) annotation(
        Line(points = {{-28, -20}, {-4, -20}}, color = {0, 79, 0}));
      connect(ball7.center, contactsBallBall2.center_a) annotation(
        Line(points = {{-68, -60}, {-48, -60}}, color = {0, 127, 0}));
      connect(contactsBallBall2.center_b, ball8.center) annotation(
        Line(points = {{-28, -60}, {-4, -60}}, color = {0, 79, 0}));
      annotation(
        experiment(StartTime = 0, StopTime = 3, Tolerance = 1e-06, Interval = 0.00015),
        Diagram(graphics = {Text(origin = {-107, 61}, extent = {{-21, 9}, {21, -9}}, textString = "no notation"), Text(origin = {-115, 22}, extent = {{-27, 12}, {27, -12}}, textString = "rotate the left ball"), Text(origin = {-132, -20}, extent = {{-44, 18}, {44, -18}}, textString = "rotate the both ball in same direction"), Line(origin = {-79.6116, 83.3043}, points = {{-27, 0}, {11, 0}}, arrow = {Arrow.None, Arrow.Filled}, arrowSize = 6), Text(origin = {-131, -58}, extent = {{-47, 16}, {47, -16}}, textString = "rotate the both ball in different direction"), Text(origin = {-85, 91}, extent = {{-45, 13}, {45, -13}}, textString = "velocity direction of all left balls"), Text(origin = {33, 91}, extent = {{-59, 15}, {59, -15}}, textString = "all right balls is stop or rotating without moving")}, coordinateSystem(extent = {{-200, -100}, {100, 100}})),
        Icon(coordinateSystem(extent = {{-200, -100}, {100, 100}})));
    end Contact_BallBall;
    
    model Contact_BallBall_BigAndSmall
      PythagoraDevice.Components.Ball ball1(D = 1, g_f = {0, 0, 0}, v_g(fixed = true), v_g_ini = {1, 0, 0}, w_g_ini = {0, 0, -3}) annotation(
        Placement(visible = true, transformation(origin = {-53, -5}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
      PythagoraDevice.Components.ContactsBallBall contactsAddRotation annotation(
        Placement(visible = true, transformation(origin = {-14, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.Ball ball2(D = 0.3, g_f = {0, 0, 0}, position_ini = {1, 0, 0}, v_g(fixed = true), v_g_ini = {0, 0, 0}) annotation(
        Placement(visible = true, transformation(origin = {20, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(ball1.center, contactsAddRotation.center_a) annotation(
        Line(points = {{-53, -5}, {-35, -5}, {-35, -6}, {-24, -6}}, color = {0, 127, 0}));
      connect(contactsAddRotation.center_b, ball2.center) annotation(
        Line(points = {{-4, -6}, {20, -6}}, color = {0, 79, 0}));
      annotation(
        experiment(StartTime = 0, StopTime = 5, Tolerance = 1e-06, Interval = 0.001));
    end Contact_BallBall_BigAndSmall;
    
    model FallingBallOntoBox "Falling the ball onto a plane"
      //  inner Modelica.Mechanics.MultiBody.World world annotation(
      //    Placement(visible = true, transformation(origin = {-80, 78}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.FixedBox fixedBoxB(height = 0.05, length = 0.2, position_ini = {0, -0.005, 0}, width = 0.01) annotation(
        Placement(visible = true, transformation(origin = {24, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.ContactsBallFixedBox contactsBallFixedBoxB annotation(
        Placement(visible = true, transformation(origin = {98, -14}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      PythagoraDevice.Components.Ball ball_B(D = 11.3e-3, m = 1.89e-3, position_ini = {0, 0.1, 0}, v_g(fixed = true), v_g_ini = {0, 0, 0}, w_g_ini = {0, 0, 0}) annotation(
        Placement(visible = true, transformation(origin = {98, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.Ball ball_A(D = 11.3e-3, color1 = {0, 255, 255}, m = 1.89e-3, position_ini = {-0.05, 0.1, 0}, v_g(fixed = true), v_g_ini = {0, 0, 0}, w_g_ini = {0, 0, -30}) annotation(
        Placement(visible = true, transformation(origin = {-44, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.ContactsBallFixedBox contactsBallFixedBoxA annotation(
        Placement(visible = true, transformation(origin = {-44, -14}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      PythagoraDevice.Components.ContactsBallBall contactsBallBall annotation(
        Placement(visible = true, transformation(origin = {28, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(ball_B.center, contactsBallFixedBoxB.center_a) annotation(
        Line(points = {{98, 28}, {98, -4}}, color = {0, 127, 0}));
      connect(contactsBallFixedBoxB.center_b, fixedBoxB.center) annotation(
        Line(points = {{98, -24}, {98, -36}, {24, -36}, {24, -50}}, color = {0, 79, 0}));
      connect(ball_A.center, contactsBallFixedBoxA.center_a) annotation(
        Line(points = {{-44, 28}, {-44, -4}}, color = {0, 127, 0}));
      connect(ball_A.center, contactsBallBall.center_a) annotation(
        Line(points = {{-44, 28}, {18, 28}, {18, 26}}, color = {0, 127, 0}));
      connect(contactsBallBall.center_b, ball_B.center) annotation(
        Line(points = {{38, 26}, {98, 26}, {98, 28}}, color = {0, 79, 0}));
      connect(contactsBallFixedBoxA.center_b, fixedBoxB.center) annotation(
        Line(points = {{-44, -24}, {-44, -36}, {24, -36}, {24, -50}}, color = {0, 79, 0}));
    protected
      annotation(
        experiment(StartTime = 0, StopTime = 3, Tolerance = 1e-06, Interval = 0.0003),
        Diagram(graphics = {Rectangle(extent = {{28, -12}, {28, -12}}), Text(origin = {-44, 49}, extent = {{-44, 5}, {44, -5}}, textString = "With rotation"), Text(origin = {98, 49}, extent = {{-44, 5}, {44, -5}}, textString = "No rotation")}, coordinateSystem(extent = {{-200, -100}, {200, 100}})),
        Icon(coordinateSystem(extent = {{-200, -100}, {200, 100}})));
    end FallingBallOntoBox;
  end Examples;

  package Components
    extends Modelica.Icons.Package;

    model Ball
      /* ---------------------------------------------
                  Basic parameters
      --------------------------------------------- */
      parameter SI.Mass m(min = 0) = 4/3*pi*(D/2)^3*rho "Mass";                                 
      parameter SI.Diameter D(min = 0) = 11.3e-3 "Diameter";
      parameter StateSelect stateSelect = StateSelect.default "Priority to use s and v as states" annotation(
        Dialog(tab = "Advanced"));
      parameter SI.Force g_f[3] = {0, -9.8, 0} "Gravity force";
      
      /* ---------------------------------------------
                  Material parameters
      --------------------------------------------- */
      parameter SI.Density rho(min = 0) = 2.5e+3 "density, glass:2.3e+3kg/m^3" annotation(
        Dialog(group = "Material"));
      parameter Real E(min = 0, unit="Pa") = 7.16e+10 "Young's modulus" annotation(
        Dialog(group = "Material"));
      parameter Real nue(min = 0) = 0.23 "Poisson's ratio" annotation(
        Dialog(group = "Material"));
    
      /* ---------------------------------------------------------------------
                  Inertia parameters
      reference https://kikyousan.com/physics/dynamics1/inertia3
      ------------------------------------------------------------------------ */
      parameter SI.Inertia I_11(min = 0) = 2 / 5 * m * (D / 2) ^ 2 "(1,1) element of inertia tensor" annotation(
        Dialog(group = "Inertia tensor (resolved in center of mass, parallel to frame_a)"));
      parameter SI.Inertia I_22(min = 0) = 2 / 5 * m * (D / 2) ^ 2 "(2,2) element of inertia tensor" annotation(
        Dialog(group = "Inertia tensor (resolved in center of mass, parallel to frame_a)"));
      parameter SI.Inertia I_33(min = 0) = 2 / 5 * m * (D / 2) ^ 2 "(3,3) element of inertia tensor" annotation(
        Dialog(group = "Inertia tensor (resolved in center of mass, parallel to frame_a)"));
      parameter SI.Inertia I_21 = 0 "(2,1) element of inertia tensor" annotation(
        Dialog(group = "Inertia tensor (resolved in center of mass, parallel to frame_a)"));
      parameter SI.Inertia I_31 = 0 "(3,1) element of inertia tensor" annotation(
        Dialog(group = "Inertia tensor (resolved in center of mass, parallel to frame_a)"));
      parameter SI.Inertia I_32 = 0 "(3,2) element of inertia tensor" annotation(
        Dialog(group = "Inertia tensor (resolved in center of mass, parallel to frame_a)"));
      parameter SI.Inertia I[3, 3] = [I_11, I_21, I_31; I_21, I_22, I_32; I_31, I_32, I_33] "inertia tensor" annotation(
        Dialog(group = "Inertia tensor (resolved in center of mass, parallel to frame_a)"));
      
      /* ---------------------------------------------
                  initial condition
      --------------------------------------------- */
      parameter SI.Angle angles_g_start[3] = {0, 0, 0} "initial angles of object in world frame" annotation(
        Dialog(group = "Initialization"));
      parameter SI.Position position_ini[3] = {0, 0, 0} "initial centor position of object in world frame" annotation(
        Dialog(group = "Initialization"));
      parameter SI.AngularVelocity w_g_ini[3] = {0, 0, 0} annotation(
        Dialog(group = "Initialization"));
      parameter SI.Velocity v_g_ini[3] = {0, 0, 0} annotation(
        Dialog(group = "Initialization"));
      
      /* ---------------------------------------------
                  Visialization
      --------------------------------------------- */
      parameter Types.Color sphereColor = Modelica.Mechanics.MultiBody.Types.Defaults.BodyColor annotation(
        Dialog(tab = "Visualization"));
      parameter String shapeType = "sphere" annotation(
          Dialog(tab = "Visualization"));
      parameter Real e_vis1[3]={D*0.003,0,0} annotation(
          Dialog(tab = "Visualization"));
      parameter Real e_vis2[3]={0,D*0.003,0} annotation(
          Dialog(tab = "Visualization"));
      parameter Real color1[3]={255, 255, 0} annotation(
          Dialog(tab = "Visualization"));
      parameter Real color2[3]={0, 0, 255} annotation(
          Dialog(tab = "Visualization"));
      parameter Real color3[3]={255, 0, 0} annotation(
          Dialog(tab = "Visualization"));
      
      /* ---------------------------------------------
                  Variables
      --------------------------------------------- */   
      //Force and moment
      SI.Force F_g[3] "Force applied to an object";
      SI.Torque M_g[3] "moment";
      
      //translational
      SI.Velocity v_g[3](start = v_g_ini) "Absolute velocity of component　in world fram";
      SI.Acceleration a_g[3](start = {0, 0, 0}) "Absolute acceleration of component　in world fram";
      SI.Position position_g[3](start = position_ini) "Absolute position of center of component　in world fram";
      
      //rotational
      SI.AngularVelocity w_g[3](start = w_g_ini) "Absolute angular velocity of frame_a resolved in world fram";
      //  SI.AngularAcceleration z_d[3] "Absolute angular acceleration of frame_a resolved in body-fixed fram";
      //  SI.Angle angles_d[3];
      SI.Angle angles_g[3](start = angles_g_start);
      Frames.Orientation R;
      
      //visualization
      Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape vis1(shapeType = shapeType, color = color1, length = D, width = D, height = D, lengthDirection = {1, 0, 0}, widthDirection = {0, 1, 0}, r_shape = r_shape+e_vis1, r = position_g, R = R);
      Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape vis2(shapeType = shapeType, color = color2, length = D, width = D, height = D, lengthDirection = {1, 0, 0}, widthDirection = {0, 1, 0}, r_shape = r_shape+e_vis2, r = position_g, R = R);
      Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape vis3(shapeType = shapeType, color = color3, length = D, width = D, height = D, lengthDirection = {1, 0, 0}, widthDirection = {0, 1, 0}, r_shape = r_shape, r = position_g, R = R);
    
      /* ---------------------------------------------
                  Port
      --------------------------------------------- */  
      PythagoraDevice.Interfaces.Center_a center annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    
    protected
      parameter Real r_shape[3] = if shapeType == "sphere" then -{1, 0, 0} * D / 2 else {0, 0, 0} "Fix origin and center of gravity of sphere";
      constant String selectShape = "ball";
      parameter PythagoraDevice.Interfaces.Geometry geometry(shape = selectShape, L = D, W = D, H = D);
      parameter PythagoraDevice.Interfaces.Material material(rho = rho, E = E, nue = nue);
    
    equation
//Force and moment applied to an object
      F_g = center.f + m * g_f;
  M_g = center.t;
      
//equation of equilibrium
      F_g = m * der(v_g);
      M_g = I * der(w_g);
      
//translational
      a_g = der(v_g);
      v_g = der(position_g);
      
//rotational
      w_g = der(angles_g);
      
//Orientation instance R
      R = Frames.axesRotations(sequence = {1, 2, 3}, angles = angles_g, der_angles = w_g);
//port
      center.position = position_g;
      center.R = R;
      center.m = m;
      center.geometry = geometry;
      center.material = material;
      annotation(
        Icon(graphics = {Ellipse(fillColor = {170, 255, 255}, fillPattern = FillPattern.Sphere, extent = {{-100, 100}, {100, -100}}, endAngle = 360), Line(origin = {179.403, 0.196862}, points = {{-72, 78}, {-72, -82}}, thickness = 1.5, arrow = {Arrow.None, Arrow.Filled}, arrowSize = 18), Text(origin = {98, -86}, extent = {{-16, 8}, {36, -22}}, textString = "g"), Text(origin = {-2, 117}, extent = {{-68, 17}, {68, -17}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)),
        __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"),
        experiment(StartTime = 0, StopTime = 3, Tolerance = 1e-06, Interval = 0.006));
    end Ball;

    model FixedBox
      parameter SI.Density rho(min = 0) = 340 "density, glass:2.3e+3kg/m^3" annotation(
        Dialog(group = "Material"));
      parameter Real E(min = 0, unit="Pa") = 7.85e+9 "Young's modulus" annotation(
        Dialog(group = "Material"));
      parameter Real nue(min = 0) = 0.5 "Poisson's ratio" annotation(
        Dialog(group = "Material"));
      parameter SI.Length length = 1;
      parameter SI.Length width = 0.01;
      parameter SI.Length height = 1;
      parameter Real m = rho * length * width * height;
      
      parameter SI.Position position_ini[3] = {0, 0, 0} "centor of box";
      parameter Real lengthDirection[3] = {1, 0, 0};
      parameter Real widthDirection[3] = {0, 1, 0};
      
      //Animation
      Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape vis(shapeType = selectShape, color = {255, 85, 127}, length = length, width = width, height = height, lengthDirection = {1, 0, 0}, widthDirection = {0, 1, 0}, r_shape = r_shape, r = position_ini, R = R);
      Frames.Orientation R;
      //port
      PythagoraDevice.Interfaces.Center_a center annotation(
        Placement(visible = true, transformation(origin = {0, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    protected
      parameter Real dummy=999;
      parameter Real r_shape[3] = {-length / 2, 0, 0} "set the centor to origin";
      constant String selectShape = "box";
      PythagoraDevice.Interfaces.Geometry geometry(shape = selectShape, L = length, W = width, H = height, m=dummy);
    equation
      R = Modelica.Mechanics.MultiBody.Frames.nullRotation();
      
    //port
      center.position = position_ini;
      center.geometry = geometry;
      center.R = R;
      center.m = m;
      center.material.rho = rho;
      center.material.E = E;
      center.material.nue = nue;
      annotation(
        uses(Modelica(version = "3.2.3")),
        Icon(coordinateSystem(initialScale = 0.1), graphics = {Rectangle(fillColor = {255, 85, 127}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-100, 20}, {100, -20}})}));
    end FixedBox;

    model Wall_rotational
      /* ---------------------------------------------
                                                parameters
                                    --------------------------------------------- */
      parameter SI.Mass m(min = 0) = 1 "Mass";
      parameter SI.Diameter D(min = 0) = 1 "Diameter";
      parameter StateSelect stateSelect = StateSelect.default "Priority to use s and v as states" annotation(
        Dialog(tab = "Advanced"));
      parameter SI.Force g_f[3] = {0, -9.8, 0} "Gravity force";
      parameter SI.Density rho(min = 0) = 2.5e+3 "density, glass:2.3e+3kg/m^3" annotation(
        Dialog(group = "Material"));
      parameter Real E(min = 0, unit="Pa") = 7.16e+10 "Young's modulus" annotation(
        Dialog(group = "Material"));
      parameter Real nue(min = 0) = 0.23 "Poisson's ratio" annotation(
        Dialog(group = "Material"));
      //inertia
      /* ---------------------------------------------------------------------
                                          reference https://kikyousan.com/physics/dynamics1/inertia3
                                          ------------------------------------------------------------------------ */
      parameter SI.Inertia I_11(min = 0) = 2 / 5 * m * (D / 2) ^ 2 "(1,1) element of inertia tensor" annotation(
        Dialog(group = "Inertia tensor (resolved in center of mass, parallel to frame_a)"));
      parameter SI.Inertia I_22(min = 0) = 2 / 5 * m * (D / 2) ^ 2 "(2,2) element of inertia tensor" annotation(
        Dialog(group = "Inertia tensor (resolved in center of mass, parallel to frame_a)"));
      parameter SI.Inertia I_33(min = 0) = 2 / 5 * m * (D / 2) ^ 2 "(3,3) element of inertia tensor" annotation(
        Dialog(group = "Inertia tensor (resolved in center of mass, parallel to frame_a)"));
      parameter SI.Inertia I_21 = 0 "(2,1) element of inertia tensor" annotation(
        Dialog(group = "Inertia tensor (resolved in center of mass, parallel to frame_a)"));
      parameter SI.Inertia I_31 = 0 "(3,1) element of inertia tensor" annotation(
        Dialog(group = "Inertia tensor (resolved in center of mass, parallel to frame_a)"));
      parameter SI.Inertia I_32 = 0 "(3,2) element of inertia tensor" annotation(
        Dialog(group = "Inertia tensor (resolved in center of mass, parallel to frame_a)"));
      parameter SI.Inertia I[3, 3] = [I_11, I_21, I_31; I_21, I_22, I_32; I_31, I_32, I_33] "inertia tensor";
      //initial condition
      parameter SI.Angle angles_g_start[3] = {0, 0, 0} "initial angles of object in world frame" annotation(
        Dialog(group = "Initialization"));
      parameter SI.Position position_ini[3] = {0, 0, 0} "initial centor position of object in world frame" annotation(
        Dialog(group = "Initialization"));
      parameter SI.AngularVelocity w_d_ini[3] = {0, 0, 0} annotation(
        Dialog(group = "Initialization"));
      parameter SI.Velocity v_g_ini[3] = {0, 0, 0} annotation(
        Dialog(group = "Initialization"));
      //Visialization
      parameter Types.Color sphereColor = Modelica.Mechanics.MultiBody.Types.Defaults.BodyColor annotation(
        Dialog(tab = "Visualization"));
      parameter String shapeType = "box" annotation(
        Dialog(tab = "Visualization"));
      //port
      PythagoraDevice.Interfaces.Center_a Center_a annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      //Force and moment
      SI.Force F_g[3] "Force applied to an object";
      SI.Force F_d[3] "Force applied to an object";
      SI.Torque M_d[3] "外力によるモーメント";
      //translational
      SI.Velocity v_g[3](start = v_g_ini) "Absolute velocity of component　in world fram";
      SI.Velocity v_d[3](start = {0, 0, 0}) "Absolute velocity of component in body-fixed fram";
      SI.Acceleration a_g[3](start = {0, 0, 0}) "Absolute acceleration of component　in world fram";
      SI.Position position_g[3](start = position_ini) "Absolute position of center of component　in world fram";
      //rotational
      SI.AngularVelocity w_d[3](start = w_d_ini) "Absolute angular velocity of frame_a resolved in body-fixed fram";
      SI.AngularVelocity w_g[3] "Absolute angular velocity of frame_a resolved in world fram";
      //  SI.AngularAcceleration z_d[3] "Absolute angular acceleration of frame_a resolved in body-fixed fram";
      //  SI.Angle angles_d[3];
      SI.Angle angles_g[3](start = angles_g_start);
      Frames.Orientation R;
      //visualization
      Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape vis(shapeType = shapeType, color = {255, 255, 255}, length = D, width = D, height = D, lengthDirection = {1, 0, 0}, widthDirection = {0, 1, 0}, r_shape = r_shape, r = position_g, R = R);
    protected
      //visualization
      parameter Real r_shape[3] = -{1, 0, 0} * D / 2;
      //collision
      constant String selectShape = "ball";
      parameter PythagoraDevice.Interfaces.Geometry geometry(shape = selectShape, L = D, W = D, H = D);
    equation
//Force and moment applied to an object
      F_g = Center_a.f + m * g_f;
      M_d = Center_a.t;
//Coordinate transformation
      F_d = Frames.resolve2(R, F_g);
      v_g = Frames.resolve1(R, v_d);
      w_g = Frames.resolve1(R, w_d);
      
//equation of equilibrium
      F_d = m * (der(v_d) + cross(w_d, v_d));
      M_d = I * der(w_d) + cross(w_d, I * w_d);
//translational
      a_g = der(v_g);
      v_g = der(position_g);
//rotational
      w_g = der(angles_g);
//Orientation instance R
      R = Frames.axesRotations(sequence = {1, 2, 3}, angles = angles_g, der_angles = w_g);
//port
      Center_a.position = position_g;
      Center_a.R = R;
      Center_a.geometry = geometry;
      Center_a.material.rho = rho;
      Center_a.material.E = E;
      Center_a.material.nue = nue;
      annotation(
        uses(Modelica(version = "3.2.3")),
        Icon(coordinateSystem(initialScale = 0.1), graphics = {Rectangle(fillColor = {170, 255, 255}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-100, 20}, {100, -20}})}));
    end Wall_rotational;
  
    model MultiBalls "Test model"
      parameter Integer n = 5 "Number of balls";
      parameter SI.Position position[n, 3] = fill(0.1, n, 3) "center position of balls";
      parameter SI.Velocity v_start[n, 3] = fill(0, n, 3) "initial velocity of balls";
      parameter SI.Diameter D[n] = fill(0.1, n) "Diameter";
      parameter SI.Mass m[n] = fill(0.1, n) "Mass";
      parameter SI.Acceleration g[n, 3] = fill({0, -9.8, 0}, n);
      PythagoraDevice.Components.Ball ball[n](D = D, g_f = g, m = m, position = position, v(fixed = true, start = v_start)) annotation(
        Placement(visible = true, transformation(origin = {-2, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.Contacts contacts[num] annotation(
        Placement(visible = true, transformation(origin = {-2, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    protected
      parameter Integer num = integer(n * (n - 1) / 2);
      parameter Integer A1[num] = PythagoraDevice.Utilities.combinationMatrix2(n);
      parameter Integer A2[num] = PythagoraDevice.Utilities.combinationMatrix3(n);
      Interfaces.Center_a Center_a[n] annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      for i in 1:num loop
        connect(ball[A1[i]].Center_a, contacts[i].Center_a);
        connect(ball[A2[i]].Center_a, contacts[i].Center_b);
      end for;
      for j in 1:n loop
        connect(ball[j].Center_a, Center_a[j]);
      end for;
      annotation(
        Icon(graphics = {Ellipse(origin = {-44, 75}, fillColor = {170, 255, 0}, fillPattern = FillPattern.Sphere, extent = {{-28, 27}, {48, -47}}, endAngle = 360), Ellipse(origin = {76, 17}, fillColor = {0, 255, 255}, fillPattern = FillPattern.Sphere, extent = {{-68, 67}, {32, -31}}, endAngle = 360), Ellipse(origin = {-78, -21}, fillColor = {85, 255, 0}, fillPattern = FillPattern.Sphere, extent = {{-22, 23}, {86, -79}}, endAngle = 360), Ellipse(origin = {36, -25}, fillColor = {255, 0, 0}, fillPattern = FillPattern.Sphere, extent = {{2, -1}, {56, -53}}, endAngle = 360)}, coordinateSystem(initialScale = 0.1)),
        experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002));
    end MultiBalls;

    model FixedMultiBalls "Test model"
      parameter Integer n = 3;
      parameter Integer num = integer(n * (n - 1) / 2);
      parameter Real position[n, 3] = fill(0.1, n, 3);
      parameter Real v_start[n, 3] = fill(0, n, 3);
      parameter Real D[n] = fill(0.1, n);
      parameter Real m[n] = fill(0.1, n);
      PythagoraDevice.Components.FixedBall ball[n](sphereColor = fill({255, 255, 255}, n), D = D, m = m, position = position, s(fixed = true, start = position), v(fixed = false)) annotation(
        Placement(visible = true, transformation(origin = {-2, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.Center_a Center_a[n] annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      for j in 1:n loop
        connect(ball[j].Center_a, Center_a[j]);
        ball[j].Center_a.xyz = position[j];
      end for;
      annotation(
        Icon(graphics = {Ellipse(origin = {-44, 75}, fillColor = {85, 0, 0}, fillPattern = FillPattern.Sphere, extent = {{-28, 27}, {48, -47}}, endAngle = 360), Ellipse(origin = {76, 17}, fillPattern = FillPattern.Sphere, extent = {{-68, 67}, {32, -31}}, endAngle = 360), Ellipse(origin = {-78, -21}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Sphere, extent = {{-22, 23}, {86, -79}}, endAngle = 360), Ellipse(origin = {36, -25}, fillColor = {0, 85, 127}, fillPattern = FillPattern.Sphere, extent = {{2, -1}, {56, -53}}, endAngle = 360)}, coordinateSystem(initialScale = 0.1)),
        experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002));
    end FixedMultiBalls;

    model ContactsBallBall
      import PythagoraDevice.Utilities.*;
      extends Interfaces.Partials.TwoPort;
      
    // parameter
      parameter StateSelect stateSelect = StateSelect.prefer "Priority to use s_rel and v_rel as states" annotation(
        HideResult = true,
        Dialog(tab = "Advanced"));
      parameter SI.Distance s_nominal[3] = fill(1e-4, 3) "Nominal value of s_rel (used for scaling)" annotation(
        Dialog(tab = "Advanced"));
      parameter Real mu = 0.4;
      parameter Real mu0 = 0.6 "coefficient of static friction";
      parameter Real e=0.6 "coefficient of restitution";
      parameter Real v_brk = 1e-6;
      parameter Real d_vis = 1e-3;  
    
      
    // Variants
      SI.TranslationalSpringConstant c_t(final min = 0) "Spring constant";
      SI.TranslationalDampingConstant d_t(final min = 0) "Damping constant";
      SI.TranslationalSpringConstant c_n(final min = 0) "Spring constant";
      SI.TranslationalDampingConstant d_n(final min = 0) "Damping constant";
      SI.Position s_rel[3](nominal = s_nominal) "Relative distance (= center_b.s - center_a.s)";
      SI.Velocity v_rel[3] "Relative velocity (= der(s_rel))";
      SI.Force f[3] "Forces between Centers";
      SI.Position s_a[3];
      SI.Position s_b[3];
      Real s_penetrate[3] "relative position vector of overlapping two objects";
      
      Real s_radius[3];
      Real torque[3] "Moments between Centers";
      SI.Force f_n[3] "normal forces between Centers";
      SI.Force f_t[3] "tangent forces between Centers";
      Real s_t[3];
      SI.Velocity v_t[3];
      SI.Velocity v_n[3];
      SI.Length D_a;
      SI.Length D_b;
      Real w_a[3];
      Real w_b[3];
      Real n_ab[3] "Normal vector of between two objects from flang.a";
      Real f_t_friction[3];
      Real f_t_penetrate[3];
      Integer contactState;
      String shapeCombi;
      SI.Length L_gap;
      SI.Length L_threshold;
      SI.Length L_penetrate "penetration";
      Modelica.SIunits.Force f_nc[3] "Spring force";
      Modelica.SIunits.Force f_nd[3] "Linear damping force";
      Modelica.SIunits.Force f_nd2[3] "Linear damping force";
    
    
      Real v_gap "relative velocity; If two objects become closer, the relative velocity is negative.
    If two objects become come away, the relative velocity is positive.";
      Real E_star;
      Real R_star;
      Real m_star;
      Real alpha;
      Real E_a;
      Real E_b;
      Real nue_a;
      Real nue_b;
      Real m_a;
      Real m_b;
      Real G_star;
      
      Real f_coulomb;
      Real f_brk;
      Real v_st;
      Real v_coulomb;
      Real f_stribeck[3];
    equation
    //nonlinear spring stiffness c_n and damping coefficient d_n for normal direction
      1/E_star = (1-nue_a^2)/E_a+(1-nue_b^2)/E_b;
      1/R_star = 1/(D_a/2)+1/(D_b/2);
      1/m_star = 1/m_a + 1/m_b;
      alpha = log(e)*sqrt(1/(log(e)^2+pi^2));
      c_n = 4/3*E_star*sqrt(R_star*L_penetrate);
      d_n = 2*sqrt((m_star * c_n)/(1+(pi/log(e))^2));
      
    //nonlinear spring stiffness c_t and damping coefficient d_t for tangent direction
      c_t = 8*G_star*sqrt(R_star * L_penetrate);
      1/G_star = 2 * (2 - nue_a) * (1 + nue_a) / E_a + 2 * (2 - nue_b) * (1 + nue_b) / E_b;
      d_t = d_n;
    
    //Relative geometric relationships and penetration depth
      s_rel = s_b - s_a;
      L_gap = length(s_rel);
      L_threshold = D_a / 2 + D_b / 2;
      L_penetrate = max(L_threshold - L_gap, 0);
      n_ab = normalize(s_rel);
      s_penetrate = L_penetrate * n_ab;
      v_rel = der(s_rel);
      v_gap = der(L_gap);
    
    //normal force
      f_n = f_nc + f_nd;
      f_nc = -c_n * s_penetrate;
      f_nd = smooth(0, noEvent(if L_penetrate <> 0 then (if v_gap < 0 and length(f_nd2) > length(f_nc) then f_nc elseif v_gap >= 0 and  length(f_nd2) > length(f_nc) then -f_nc else f_nd2 ) else zeros(3)));  
      f_nd2 = -d_n*length(v_rel)*n_ab;
    
    //tangent force
      v_n = v_rel * n_ab * n_ab;
      v_t = v_rel - v_n + cross(D_a / 2 * w_a + D_b / 2 * w_b, n_ab);
      der(s_t) = v_t;
    
    /* The reference of this friciton model is https://jp.mathworks.com/help/physmod/simscape/ref/translationalfriction.html */
      v_st = v_brk * sqrt(2);
      v_coulomb = v_brk / 10;
      f_coulomb = mu * length(f_n);
      f_brk = mu0 * length(f_n);
      for i in 1:3 loop   // i shall not be declared
            f_stribeck[i] = -sqrt(2 * exp(1)) * (f_brk - f_coulomb) * exp(-(v_t[i] / v_st) ^ 2) * v_t[i] / v_st;
      end for;
      f_t_friction = f_stribeck - f_coulomb * tanh(v_t / v_coulomb);
    
      f_t_penetrate = (-c_t * s_t) - d_t * v_t;
      if L_penetrate == 0 then
        f_t = zeros(3);
        contactState = 0;
      else
        if length(f_t_penetrate) < length(f_t_friction) then
          f_t = f_t_penetrate;
          contactState = 1;
        else
          f_t = f_t_friction;
          contactState = 2;
        end if;
      end if;
    
    //contact force
      f = f_n + f_t;
      
    //contact torque
      s_radius = center_a.geometry.L * n_ab;
      torque = cross(s_radius, f_t);
    
    //port
      center_a.position = s_a;
      center_b.position = s_b;
      center_a.f = -f;
      center_b.f = f;
      center_a.t = -torque;
      center_b.t = -torque;
      center_a.R.w = w_a;
      center_b.R.w = w_b;
      center_a.geometry.L = D_a;
      center_b.geometry.L = D_b;
      center_a.material.E = E_a;
      center_b.material.E = E_b;
      center_a.material.nue = nue_a;
      center_b.material.nue = nue_b;
      center_a.m = m_a;
      center_b.m = m_b;  
    
    //dummy
      shapeCombi = PythagoraDevice.Components.Functions.determineShapeCombination(center_a.geometry.shape, center_b.geometry.shape);
      annotation(
        Icon(graphics = {Polygon(origin = {-1, 9}, fillColor = {255, 255, 0}, fillPattern = FillPattern.Solid, points = {{-55, 87}, {-49, 43}, {-47, 21}, {-95, 15}, {-47, -29}, {-59, -87}, {-13, -37}, {-11, -79}, {5, -47}, {21, -89}, {39, -37}, {87, -65}, {51, -15}, {93, -1}, {53, 15}, {77, 27}, {61, 27}, {85, 79}, {23, 43}, {9, 89}, {-7, 45}, {-55, 87}})}),
        experiment(StartTime = 0, StopTime = 0.74, Tolerance = 1e-6, Interval = 1.48e-05));
    end ContactsBallBall;
    
    model ContactsBallFixedBox
      extends PythagoraDevice.Interfaces.Partials.PartialContact;
      
      /* ---------------------------------------------
                  Parameters
      --------------------------------------------- */
      parameter Real eps(min=0) = 1e-12 "minimum value of L_gap. To avoid division by zero where v_gap = der(L_gap), eps introduced." annotation(
        Dialog(group = "Advanced"));
    
      /* ---------------------------------------------
                  Force and torque variables
      --------------------------------------------- */   
      SI.Force f_coulomb "Coulomb friction";
      SI.Force f_brk "breakaway friction";
      SI.Force f_stribeck[3] "Stribeck friction";
        
      /* ---------------------------------------------------------------------------
                  Position, velocity, length and angular velocity variables
      ---------------------------------------------------------------------------- */  
      SI.Velocity v_st "Stribeck velocity threshold";
      SI.Velocity v_coulomb "Coulomb velocity threshold";
      SI.Velocity v_a[3] "velocity of conneting port a ";
      SI.Position s_closestPoint[3] "closest point to ball on box in global frame";
      SI.Position XYZmin[3];
      SI.Position XYZmax[3];
      SI.Length D_a;
      
    equation
    //nonlinear spring stiffness c_n and damping coefficient d_n for normal direction
      R_star = D_a / 2;
      m_star = m_a;
      
    //Relative geometric relationships and penetration depth
      //closest point of two objects
      XYZmin = s_b - {center_b.geometry.L, center_b.geometry.W, center_b.geometry.H} / 2;
      XYZmax = s_b + {center_b.geometry.L, center_b.geometry.W, center_b.geometry.H} / 2;
      for i in 1:3 loop
        s_closestPoint[i] = max(XYZmin[i], min(s_a[i], XYZmax[i]));
      end for;
      
      //Relative geometric relationships
      s_rel = s_closestPoint - s_a;
      L_gap = max(length(s_rel),eps); 
      L_threshold = D_a / 2;
      v_rel = der(s_a);
      n_ab = normalize(s_rel);
      
    //tangent force
    /* The reference of this friciton model is https://jp.mathworks.com/help/physmod/simscape/ref/translationalfriction.html */
      v_st = v_brk * sqrt(2);
      v_coulomb = v_brk / 10;
      f_coulomb = mu * length(f_n);
      f_brk = mu0 * length(f_n);
      for i in 1:3 loop   // i shall not be declared
            f_stribeck[i] = -sqrt(2 * exp(1)) * (f_brk - f_coulomb) * exp(-(v_t[i] / v_st) ^ 2) * v_t[i] / v_st;
      end for;
      f_t_friction = f_stribeck - f_coulomb * tanh(v_t / v_coulomb);
      L_contact_a = D_a / 2;
      L_contact_b = length(s_closestPoint - s_b);
    
    //contact torque
      s_torque = D_a * n_ab;
     
    //viscous friciton
      v_a = der(s_a);
      if L_penetrate == 0 then
        f_vis = zeros(3);
      else
        f_vis = -d_vis * v_a;
      end if;
      
    //port
      center_a.geometry.L = D_a;
    
    //assertion
      assert(f_brk >= f_coulomb, "f_brk must be greater than f_coulomb.");
    end ContactsBallFixedBox;

    package Functions
      function determineShapeCombination
        // shapeCombination
        input String shape_a;
        input String shape_b;
        output String shapeCombi;
      algorithm
        if shape_a == "ball" and shape_b == "ball" then
          shapeCombi := "ball-ball";
        elseif shape_a == "ball" and shape_b == "box" or shape_a == "box" and shape_b == "ball" then
          shapeCombi := "ball-box";
        else
          shapeCombi := "other";
        end if;
      end determineShapeCombination;
    end Functions;
    annotation(
      Icon(graphics = {Ellipse(origin = {-46, 18}, fillColor = {0, 255, 0}, fillPattern = FillPattern.Sphere, extent = {{-50, 52}, {50, -52}}, endAngle = 360), Rectangle(origin = {64, -20}, fillColor = {170, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-22, 74}, {22, -74}})}));
  end Components;

  package Utilities
    extends Modelica.Icons.UtilitiesPackage;

    function cosTwoVectors
      input Real v1[:];
      input Real v2[:];
      output Real cos;
    algorithm
      cos := v1 * v2 / (Modelica.Math.Vectors.length(v1) * Modelica.Math.Vectors.length(v2));
    end cosTwoVectors;

    function combinationMatrix2
      input Integer n;
      output Integer A1[num];
      parameter Integer num = integer(n * (n - 1) / 2);
      Integer f = n;
      Integer k;
      Integer m;
    algorithm
      k := 0;
      for i in 1:n - 1 loop
        f := f - 1;
        m := i;
        for j in f:(-1):1 loop
          k := k + 1;
          A1[k] := i;
        end for;
      end for;
    end combinationMatrix2;

    function combinationMatrix3
      input Integer n;
      output Integer A2[num];
      parameter Integer num = integer(n * (n - 1) / 2);
      Integer f = n;
      Integer k;
      Integer m;
    algorithm
      k := 0;
      for i in 1:n - 1 loop
        f := f - 1;
        m := i;
        for j in f:(-1):1 loop
          k := k + 1;
          m := m + 1;
          A2[k] := m;
        end for;
      end for;
    end combinationMatrix3;
  end Utilities;

  package Interfaces "Connectors and partial models"
    extends Modelica.Icons.InterfacesPackage;

    connector Center "One-dimensional translational Center"
      SI.Position position[3] "Absolute position of center of object";
      Frames.Orientation R "Orientation object to rotate the world frame into the body frame";
      flow SI.Force f[3] "Force resolved in world frame";
      flow SI.Torque t[3] "Torque resolved in body frame";
      Geometry geometry "Geometry infomation of object";
      Material material "Material infomation of object";
      SI.Mass m;
      annotation(
        Documentation(info = "<html>
    <p>
    This is a connector for 1D translational mechanical systems.
    It has no icon definition and is only used by inheritance from
    Center connectors to define different icons.
    </p>
    <p>
    The following variables are defined in this connector:
    </p>
    
    <blockquote><pre>
    s: Absolute position of the Center in [m]. A positive translation
    means that the Center is translated along the Center axis.
    f: Cut-force in direction of the Center axis in [N].
    </pre></blockquote>
    </html>"));
    end Center;

    connector Center_a "One-dimensional translational Center (left, Center axis directed INTO cut plane)"
      extends Center;
      annotation(
        defaultComponentName = "Center_a",
        Icon(coordinateSystem(initialScale = 0.1), graphics = {Rectangle(lineColor = {0, 127, 0}, fillColor = {0, 97, 0}, fillPattern = FillPattern.Solid, extent = {{-100, -100}, {100, 100}})}),
        Diagram(coordinateSystem(initialScale = 0.1), graphics = {Rectangle(lineColor = {0, 127, 0}, fillColor = {0, 127, 0}, fillPattern = FillPattern.Solid, extent = {{-40, -40}, {40, 40}}), Text(lineColor = {0, 127, 0}, extent = {{-160, 110}, {40, 50}}, textString = "%name")}));
    end Center_a;

    connector Center_b "One-dimensional translational Center (right, Center axis directed OUT OF cut plane)"
      extends Center;
      annotation(
        defaultComponentName = "Center_b",
        Icon(coordinateSystem(initialScale = 0.1), graphics = {Rectangle(lineColor = {0, 79, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, -100}, {100, 100}})}),
        Diagram(coordinateSystem(initialScale = 0.1), graphics = {Rectangle(lineColor = {0, 127, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-40, -40}, {40, 40}}), Text(lineColor = {0, 127, 0}, extent = {{-40, 110}, {160, 50}}, textString = "%name")}));
    end Center_b;

    record Material
      SI.Density rho(min = 0) "density";
      Real E(min = 0, unit="Pa") "Young's modulus";
      Real nue(min = 0) "Poisson's ratio";
    end Material;

    record Geometry
      String shape "shape of object. ex. ball, box";
      SI.Length L "length";
      SI.Length W "wide";
      SI.Length H "height";
      SI.Mass m "Mass";  
    end Geometry;

    package Partials "Partial models"
      extends Modelica.Icons.BasesPackage;

      partial model TwoPort
        //port
        PythagoraDevice.Interfaces.Center_a center_a annotation(
          Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        PythagoraDevice.Interfaces.Center_b center_b annotation(
          Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation

      end TwoPort;

      partial model PartialContact
      
        import PythagoraDevice.Utilities.*;
        extends PythagoraDevice.Interfaces.Partials.TwoPort;
        
        /* ---------------------------------------------
                    Parameters
        --------------------------------------------- */
        parameter SI.Distance s_nominal[3] = fill(1e-4, 3) "Nominal value of s_rel (used for scaling)" annotation(
          Dialog(tab = "Advanced"));
        parameter Real mu = 0.4 "Coefficient of dynamic friction";
        parameter Real mu0 = 0.6 "Coefficient of static friction";
        parameter Real e = 0.6 "Coefficient of restitution";
        parameter SI.Velocity v_brk(min=0) = 1e-6 "Breakaway friction velocity";
        parameter SI.TranslationalDampingConstant d_vis(min=0) = 1e-3 "Viscous friction coefficient";
      
        /* ---------------------------------------------
                    Force and torque variables
        --------------------------------------------- */    
        SI.Force f[3] "Forces between Centers";
        SI.Force f_n[3] "Forces between Centers";
        SI.Force f_t[3] "Forces between Centers";
        SI.Force f_nc[3] "Spring force";
        SI.Force f_nd[3] "Linear damping force";
        SI.Force f_nd2[3] "Linear damping force: dummy";
        SI.Force f_t_friction[3] "Friction force";
        SI.Force f_t_penetrate[3] "Tangent force by penetration";
        SI.Force f_vis[3] "viscosity force";
        SI.Torque torque[3] "Moments between Centers";
      
        /* ---------------------------------------------------------------------------
                    Position, velocity, length and angular velocity variables
        ---------------------------------------------------------------------------- */    
        SI.Position s_rel[3](nominal = s_nominal) "Relative distance (= center_b.s - center_a.s)";
        SI.Position s_penetrate[3] "Relative position vector of overlapping two objects";
        SI.Position s_torque[3] "Vector of between ball(center_a port) to contact point";
        SI.Position s_a[3] "Position of center_a port";
        SI.Position s_b[3] "Position of center_b port";
        SI.Position s_t[3] "tangent direction displacement";
        SI.Velocity v_rel[3] "Relative velocity (= der(s_rel))";
        SI.Velocity v_n[3] "normal direction velocity";
        SI.Velocity v_t[3] "tangent direction velocity";
        SI.Velocity v_gap "relative velocity; If two objects become closer, the relative velocity is negative.
              If two objects become come away, the relative velocity is positive.";
        SI.Length L_gap "gap of two objects";
        SI.Length L_threshold "Collision detection threshold";
        SI.Length L_penetrate "penetration";
        SI.Length L_contact_a "length between center of gravity to contact point on the object A";
        SI.Length L_contact_b "length between center of gravity to contact point on the object B";
        SI.AngularVelocity w_a[3];
        SI.AngularVelocity w_b[3];
        Real n_ab[3] "Normal vector of between two objects from flang.a";
        
        /* ---------------------------------------------
                    Spring and damper variables
                    Reference: http://penguinitis.g1.xrea.com/study/note/discrete_element_method.pdf
        --------------------------------------------- */  
        SI.TranslationalSpringConstant c_n(final min = 0) "Spring constant";
        SI.TranslationalSpringConstant c_t(final min = 0) "Spring constant";
        SI.TranslationalDampingConstant d_n(final min = 0) "Damping constant";
        SI.TranslationalDampingConstant d_t(final min = 0) "Damping constant";
        Real E_star;
        Real R_star;
        Real m_star;
        Real alpha;
        Real E_a;
        Real E_b;
        Real nue_a;
        Real nue_b;
        Real m_a;
        Real m_b;
        Real G_star;
      
      //Comfirm state
        Integer contactState;
        String shapeCombi;
        
      equation
      /* ---------------------------------------------
      Notice: Formulas that need to be defined at the inheritance destination are commented out
      --------------------------------------------- */
      
      //nonlinear spring stiffness c_n and damping coefficient d_n for normal direction
        1/E_star = (1-nue_a^2)/E_a+(1-nue_b^2)/E_b;
      //  1/R_star = 1/(D_a/2)+1/(D_b/2);
      //  1/m_star = 1/m_a + 1/m_b;
        alpha = log(e)*sqrt(1/(log(e)^2+pi^2));
        c_n = 4/3*E_star*sqrt(R_star*L_penetrate);
        d_n = 2*sqrt((m_star * c_n)/(1+(pi/log(e))^2));
        
      //nonlinear spring stiffness c_t and damping coefficient d_t for tangent direction
        c_t = 8*G_star*sqrt(R_star * L_penetrate);
        1/G_star = 2 * (2 - nue_a) * (1 + nue_a) / E_a + 2 * (2 - nue_b) * (1 + nue_b) / E_b;
        d_t = d_n;
      
      //Relative geometric relationships and penetration depth
      //  s_rel = s_b - s_a;
      //  L_gap = length(s_rel);
      //  L_threshold = D_a / 2 + D_b / 2;
      //  n_ab = normalize(s_rel);
      //  v_rel = der(s_rel);
        L_penetrate = max(L_threshold - L_gap, 0);
        s_penetrate = L_penetrate * n_ab;
        v_gap = der(L_gap);
      
      //normal force
        f_n = f_nc + f_nd;
        f_nc = -c_n * s_penetrate;
        f_nd = smooth(0, noEvent(if L_penetrate <> 0 then (if v_gap < 0 and length(f_nd2) > length(f_nc) then f_nc elseif v_gap >= 0 and  length(f_nd2) > length(f_nc) then -f_nc else f_nd2 ) else zeros(3)));  
        f_nd2 = -d_n*v_rel*n_ab*n_ab;
      
      //tangent force
        v_n = v_rel * n_ab * n_ab;
        v_t = v_rel - v_n + cross(L_contact_a * w_a + L_contact_b * w_b, n_ab);
        der(s_t) = v_t;
      
      //  f_t_friction = -mu * length(f_n) * normalize(v_t);
        f_t_penetrate = (-c_t * s_t) - d_t * v_t;
        if L_penetrate == 0 then
          f_t = zeros(3);
          contactState = 0;
        else
          if length(f_t_penetrate) < length(f_t_friction) then
            f_t = f_t_penetrate;
            contactState = 1;
          else
            f_t = f_t_friction;
            contactState = 2;
          end if;
        end if;
      
      //contact force
        f = f_n + f_t + f_vis;
        
      //contact torque
      //  s_torque = center_a.geometry.L * n_ab;
        torque = cross(s_torque, f_t);
      
      //port
        center_a.position = s_a;
        center_b.position = s_b;
        center_a.f = -f;
        center_b.f = f;
        center_a.t = -torque;
        center_b.t = -torque;
        center_a.R.w = w_a;
        center_b.R.w = w_b;
        center_a.material.E = E_a;
        center_b.material.E = E_b;
        center_a.material.nue = nue_a;
        center_b.material.nue = nue_b;
        center_a.m = m_a;
        center_b.m = m_b;  
      
      //dummy
        shapeCombi = PythagoraDevice.Components.Functions.determineShapeCombination(center_a.geometry.shape, center_b.geometry.shape);
        annotation(
          Icon(graphics = {Polygon(origin = {-1, 9}, fillColor = {255, 255, 0}, fillPattern = FillPattern.Solid, points = {{-55, 87}, {-49, 43}, {-47, 21}, {-95, 15}, {-47, -29}, {-59, -87}, {-13, -37}, {-11, -79}, {5, -47}, {21, -89}, {39, -37}, {87, -65}, {51, -15}, {93, -1}, {53, 15}, {77, 27}, {61, 27}, {85, 79}, {23, 43}, {9, 89}, {-7, 45}, {-55, 87}})}),
          experiment(StartTime = 0, StopTime = 0.74, Tolerance = 1e-6, Interval = 1.48e-05));
      end PartialContact;
    end Partials;
  end Interfaces;

  model Enviroment //Test model
    parameter SI.Acceleration g[3] = {0, -Modelica.Constants.g_n, 0} "Constant gravity acceleration" annotation(
      Dialog(group = "Gravity"));
    parameter Real wind;
    parameter Real[3] windDirection = {1, 0, 0};
  equation

    annotation(
      defaultComponentName = "enviroment",
      defaultComponentPrefixes = "inner",
      missingInnerMessage = "No \"world\" component is defined. A default world
  component with the default gravity field will be used
  (g=9.81 in negative y-axis). If this is not desired,
  drag Modelica.Mechanics.MultiBody.World into the top level of your model.",
      Icon(graphics = {Text(origin = {-17, 99}, extent = {{-43, 17}, {83, -33}}, textString = "g"), Text(origin = {-2, -138}, extent = {{-96, 26}, {94, -8}}, textString = "%name"), Rectangle(origin = {0, 20}, fillColor = {0, 170, 255}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{-10, 40}, {12, -50}}), Polygon(origin = {-20, -59}, fillColor = {0, 170, 255}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, points = {{-20, 39}, {20, -19}, {60, 39}, {-20, 39}}), Line(origin = {-179.87, -99.53}, points = {{278.5, 0}, {80.5, 0}, {80.5, 0}}, thickness = 3, arrow = {Arrow.Open, Arrow.None}, arrowSize = 10), Line(origin = {-292.487, -63.202}, points = {{194.5, 162}, {194.5, -38}, {194.5, -38}}, thickness = 3, arrow = {Arrow.Open, Arrow.None}, arrowSize = 10), Text(origin = {-85, 106}, extent = {{-25, 18}, {55, -44}}, textString = "y"), Text(origin = {75, -52}, extent = {{-25, 18}, {55, -44}}, textString = "x")}, coordinateSystem(initialScale = 0.1)));
  end Enviroment;

  package Sources
    extends Modelica.Icons.SourcesPackage;

    model InputForceTorque
      PythagoraDevice.Interfaces.Center_a Center_a annotation(
        Placement(visible = true, transformation(origin = {86, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      parameter SI.Force Fx = 0 "external x direction force in world flame";
      parameter SI.Force Fy = 0 "external y direction force in world flame";
      parameter SI.Force Fz = 0 "external z direction force in world flame";
      parameter SI.Torque tx = 0 "external x direction torque in body flame";
      parameter SI.Torque ty = 0 "external y direction torque in body flame";
      parameter SI.Torque tz = 0 "external z direction torque in body flame";
    equation
      Center_a.f = -{Fx, Fy, Fz};
      Center_a.t = -{tx, ty, tz};
      annotation(
        Icon(graphics = {Rectangle(origin = {0, 30}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 30}, {100, -30}}), Rectangle(origin = {0, -30}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 30}, {100, -30}}), Text(origin = {-15.938, 37.5082}, extent = {{-48.0619, 22.4918}, {77.9386, -33.5082}}, textString = "Force"), Text(origin = {-16.784, -18.918}, extent = {{-57.2165, 28.918}, {92.7835, -43.082}}, textString = "Torque")}, coordinateSystem(extent = {{-100, -60}, {100, 60}})),
        Diagram(coordinateSystem(extent = {{-100, -60}, {100, 60}})));
    end InputForceTorque;

    model InputForceTorqueExternal
      PythagoraDevice.Interfaces.Center_a Center_a annotation(
        Placement(visible = true, transformation(origin = {86, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput t[3] annotation(
        Placement(visible = true, transformation(origin = {-120, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-120, -38}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput f[3] annotation(
        Placement(visible = true, transformation(origin = {-120, 80}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-120, 40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    equation
      Center_a.f = -f;
      Center_a.t = -t;
      annotation(
        Icon(graphics = {Rectangle(origin = {0, 30}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 30}, {100, -30}}), Rectangle(origin = {0, -30}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 30}, {100, -30}}), Text(origin = {-15.938, 37.5082}, extent = {{-48.0619, 22.4918}, {77.9386, -33.5082}}, textString = "Force"), Text(origin = {-16.784, -18.918}, extent = {{-57.2165, 28.918}, {92.7835, -43.082}}, textString = "Torque")}, coordinateSystem(extent = {{-100, -60}, {100, 60}})),
        Diagram(coordinateSystem(extent = {{-100, -60}, {100, 60}})));
    end InputForceTorqueExternal;
  end Sources;
  annotation(
    uses(Modelica(version = "3.2.3")),
    Icon(graphics = {Ellipse(fillColor = {165, 244, 121}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}, endAngle = 360), Text(origin = {4, -12}, fillColor = {255, 255, 0}, fillPattern = FillPattern.Solid, extent = {{-83, 71}, {81, -69}}, textString = "ピ", fontName = "UD デジタル 教科書体 NK-B")}, coordinateSystem(initialScale = 0.1)));
end PythagoraDevice;
