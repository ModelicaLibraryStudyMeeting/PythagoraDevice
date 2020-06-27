package PythagoraDevice
  import SI = Modelica.SIunits;
  import Modelica.Math.Vectors.*;

  package Visualizers
    extends Modelica.Icons.Package;

    model Sphere
      import Modelica.Mechanics.MultiBody.Types;
      import Modelica.Mechanics.MultiBody.Interfaces;
      import SI = Modelica.SIunits;
      import Modelica.Mechanics.MultiBody.Frames;
      Frames.Orientation R;
      parameter Types.Color sphereColor = {255, 255, 204};
      parameter SI.Diameter sphereDiameter = 0.1;
      Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape sphere(shapeType = "sphere", color = sphereColor, length = sphereDiameter, width = sphereDiameter, height = sphereDiameter, lengthDirection = {1, 0, 0}, widthDirection = {0, 1, 0}, r_shape = -{1, 0, 0} * sphereDiameter / 2, r = r_0, R = R);
      parameter StateSelect stateSelect = StateSelect.avoid;
      SI.Position r_0[3](start = {0, 0, 0}, each stateSelect = stateSelect);
      Modelica.Blocks.Interfaces.RealInput u[3] annotation(
        Placement(visible = true, transformation(origin = {-114, 6}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-114, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    equation
      r_0 = u;
      R.T = identity(3);
      R.w = zeros(3);
      annotation(
        uses(Modelica(version = "3.2.3")),
        Icon(graphics = {Ellipse(fillColor = {0, 170, 255}, fillPattern = FillPattern.Sphere, extent = {{-100, 100}, {100, -100}}, endAngle = 360)}));
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

    model ContactBall
      PythagoraDevice.Components.Ball topBall(g_f = {0, 0, 0}, sphereColor = {255, 255, 255}) annotation(
        Placement(visible = true, transformation(origin = {0, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.Ball bottomBall(g_f = {0, 0, 0}, position = {0, -1.5, 0}, s(each fixed = true), v(each fixed = true, start = {0, 1, 0})) annotation(
        Placement(visible = true, transformation(origin = {0, -66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.Contacts contacts(c = 10000, d = 1) annotation(
        Placement(visible = true, transformation(origin = {0, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      inner Modelica.Mechanics.MultiBody.World world annotation(
        Placement(visible = true, transformation(origin = {-86, -24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(bottomBall.flange_a, contacts.flange_a) annotation(
        Line(points = {{0, -66}, {0, -66}, {0, -38}, {0, -38}}, color = {0, 127, 0}));
      connect(topBall.flange_a, contacts.flange_b) annotation(
        Line(points = {{0, 2}, {0, 2}, {0, -18}, {0, -18}}, color = {0, 127, 0}));
      annotation(
        experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002));
    end ContactBall;

    model FixedWallAndBall
      inner Modelica.Mechanics.MultiBody.World world annotation(
        Placement(visible = true, transformation(origin = {-86, -24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.FixedaWall fixedaWall(height = 5, length = 5, width = 0.1) annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.Ball ball(D = 1, g_f = {0, -9.8, 0}, m = 0.01, position = {0, 2, 0}, v(fixed = true, start = {1.05, 0, 0})) annotation(
        Placement(visible = true, transformation(origin = {0, 56}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.Contacts contacts(c = 100, d = 0.02) annotation(
        Placement(visible = true, transformation(origin = {0, 24}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    equation
      connect(contacts.flange_a, fixedaWall.flange_b) annotation(
        Line(points = {{0, 14}, {0, 14}, {0, 0}, {0, 0}, {0, 0}}, color = {0, 127, 0}));
      connect(contacts.flange_b, ball.flange_a) annotation(
        Line(points = {{0, 34}, {0, 34}, {0, 56}, {0, 56}}, color = {0, 79, 0}));
    protected
      annotation(
        experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-06, Interval = 0.002));
    end FixedWallAndBall;

    model ContactBall2
      PythagoraDevice.Components.Ball topBall(D = 1.5, g_f = {0, -9.8, 0}, m = 10, sphereColor = {255, 255, 255}) annotation(
        Placement(visible = true, transformation(origin = {0, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.Ball bottomBall(D = 0.6, g_f = {0, -9.8, 0}, position = {0.8, -1.3, 0}, s(each fixed = true), v(each fixed = true, start = {0, 10, 0})) annotation(
        Placement(visible = true, transformation(origin = {0, -66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.Contacts contacts(c = 1000000000, d = 1) annotation(
        Placement(visible = true, transformation(origin = {0, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      inner Modelica.Mechanics.MultiBody.World world annotation(
        Placement(visible = true, transformation(origin = {-86, -24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(bottomBall.flange_a, contacts.flange_a) annotation(
        Line(points = {{0, -66}, {0, -66}, {0, -38}, {0, -38}}, color = {0, 127, 0}));
      connect(topBall.flange_a, contacts.flange_b) annotation(
        Line(points = {{0, 2}, {0, 2}, {0, -18}, {0, -18}}, color = {0, 127, 0}));
      annotation(
        experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002));
    end ContactBall2;

    model BallInBox
      parameter Integer n = 2 "number of moving Ball";
      parameter Integer n_fix = 4 "number of fix ball";
      parameter Integer m = 3 "number of Wall";
      PythagoraDevice.Components.MultiBalls moveBalls(D = fill(0.1, n), m = fill(1, n), n = n, position = {{-0.4, 0.9, 0}, {-0.3, 0.9, 0}}, v_start = fill(0, n, 3)) annotation(
        Placement(visible = true, transformation(origin = {-2, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
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
          connect(moveBalls.flange_a[i], contacts[(i - 1) * m + j].flange_b);
          connect(fixedaWall[j].flange_b, contacts[(i - 1) * m + j].flange_a);
        end for;
      end for;
//fixBall-movingBall
      for k in 1:n loop
        for p in 1:n_fix loop
          connect(moveBalls.flange_a[k], contactsFix[(k - 1) * n_fix + p].flange_b);
          connect(fixBalls.flange_a[p], contactsFix[(k - 1) * n_fix + p].flange_a);
        end for;
      end for;
      annotation(
        experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-06, Interval = 0.002),
        __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian -d=bltdump -d=bltdump -d=bltdump ",
        __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "rungekutta"));
    end BallInBox;
  end Examples;

  package Components
    extends Modelica.Icons.Package;

    model FixedaWall
      parameter SI.Length length = 1;
      parameter SI.Length width = 0.01;
      parameter SI.Length height = 1;
      parameter SI.Position s0[3] = {0, 0, 0} "centor of box";
      //Animation
      Modelica.Blocks.Interfaces.RealOutput r[3];
      parameter Real lengthDirection[3] = {1, 0, 0};
      parameter Real widthDirection[3] = {0, 1, 0};
      PythagoraDevice.Interfaces.Flange_b flange_b annotation(
        Placement(visible = true, transformation(origin = {0, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Visualizers.Box box1(Length = length, height = height, width = width, lengthDirection = lengthDirection, widthDirection = widthDirection, r_shape = r_shape, R = Modelica.Mechanics.MultiBody.Frames.nullRotation()) annotation(
        Placement(visible = true, transformation(origin = {62, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SI.Position centerXYZ[3];
    protected
      parameter Real r_shape[3] = {-length / 2, 0, 0} "set the centor to origin";
      constant String selectShape = "box";
      PythagoraDevice.Interfaces.Geometry geometry(shape = selectShape, L = length, W = width, H = height);
    equation
      r = s0;
      connect(box1.r_0, r);
      centerXYZ = s0 + r_shape;
//port
      flange_b.xyz = s0;
      flange_b.geometry = geometry;
      annotation(
        uses(Modelica(version = "3.2.3")),
        Icon(coordinateSystem(initialScale = 0.1), graphics = {Rectangle(fillColor = {255, 85, 127}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-100, 20}, {100, -20}})}));
    end FixedaWall;

    model Ball
      import Modelica.Mechanics.MultiBody.Types;
      //parameter
      parameter SI.Mass m(min = 0) = 1 "Mass";
      parameter SI.Diameter D(min = 0) = 1 "Diameter";
      parameter Real[3] position = {0, 0, 0} "Centor position of object";
      parameter StateSelect stateSelect = StateSelect.default "Priority to use s and v as states" annotation(
        Dialog(tab = "Advanced"));
      parameter Types.Color sphereColor = Modelica.Mechanics.MultiBody.Types.Defaults.BodyColor;
      parameter SI.Force g_f[3] = {0, -9.8, 0} "Gravity force";
      //physical variables
      SI.Velocity v[3](start = {0, 0, 0}) "Absolute velocity of component";
      SI.Acceleration a[3](start = {0, 0, 0}) "Absolute acceleration of component";
      SI.Position s[3](start = {0, 0, 0}) "Absolute position of center of component";
      //Animation variables
      Modelica.Blocks.Interfaces.RealOutput xyz[3] "Postion of ball";
      PythagoraDevice.Visualizers.Sphere sphere1(sphereColor = sphereColor, sphereDiameter = D) annotation(
        Placement(visible = true, transformation(origin = {14, -54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      //port
      PythagoraDevice.Interfaces.Flange_a flange_a annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    protected
      constant String selectShape = "ball";
      PythagoraDevice.Interfaces.Geometry geometry(shape = selectShape, L = D, W = D, H = D);
    equation
      v = der(s);
      a = der(v);
      m * a = flange_a.f + m * g_f;
      xyz = s + position;
      connect(sphere1.u, xyz);
//port
      flange_a.xyz = xyz;
      flange_a.geometry = geometry;
      annotation(
        Icon(graphics = {Ellipse(fillColor = {170, 255, 255}, fillPattern = FillPattern.Sphere, extent = {{-100, 100}, {100, -100}}, endAngle = 360), Line(origin = {179.403, 0.196862}, points = {{-72, 78}, {-72, -82}}, thickness = 1.5, arrow = {Arrow.None, Arrow.Filled}, arrowSize = 18), Text(origin = {98, -86}, extent = {{-16, 8}, {36, -22}}, textString = "g"), Text(origin = {-2, 117}, extent = {{-68, 17}, {68, -17}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)),
        __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"),
        experiment(StartTime = 0, StopTime = 3, Tolerance = 1e-6, Interval = 0.006));
    end Ball;

    model Contacts
      import PythagoraDevice.Utilities.*;
      extends Interfaces.Partials.TwoPort;
      parameter SI.TranslationalSpringConstant c(final min = 0) = 100000 "Spring constant";
      parameter SI.TranslationalDampingConstant d(final min = 0) = 10 "Damping constant";
      parameter StateSelect stateSelect = StateSelect.prefer "Priority to use s_rel and v_rel as states" annotation(
        HideResult = true,
        Dialog(tab = "Advanced"));
      parameter SI.Distance s_nominal[3] = fill(1e-4, 3) "Nominal value of s_rel (used for scaling)" annotation(
        Dialog(tab = "Advanced"));
      SI.Position s_rel[3](nominal = s_nominal) "Relative distance (= flange_b.s - flange_a.s)";
      SI.Velocity v_rel[3] "Relative velocity (= der(s_rel))";
      SI.Force f[3] "Forces between flanges";
      Real s_a[3];
      Real s_b[3];
      Real vec_p[3] "penetration vector";
      Boolean contact;
      //protected
      Modelica.SIunits.Force f_c[3] "Spring force";
      Modelica.SIunits.Force f_d2[3] "Linear damping force";
      Modelica.SIunits.Force f_d[3] "Linear damping force which is limited by spring force (|f_d| <= |f_c|)";

      function shapeCombination
      /******************************************************************
      shapeCombination
      ******************************************************************/
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
      end shapeCombination;

      function penetration "return penetration[3] between two objects"
      /******************************************************************
      penetration
      ******************************************************************/
        input Real s_a[3];
        input Real s_b[3];
        input PythagoraDevice.Interfaces.Geometry geo_a;
        input PythagoraDevice.Interfaces.Geometry geo_b;
        output Real vec_p[3];
        output Boolean contact;
      protected
        Real s_rel[3];
        Real v_rel[3];
        String shape_a;
        String shape_b;
        Real gap;
        Real contactThreshold;
        Real L_p;
        Real vec_c[3];
        String shapeCombi;
      algorithm
        shape_a := geo_a.shape;
        shape_b := geo_b.shape;
        shapeCombi := shapeCombination(shape_a, shape_b);
        if shapeCombi == "ball-ball" then
          s_rel := s_a - s_b;
          gap := length(s_a - s_b);
          contactThreshold := geo_a.L / 2 + geo_b.L / 2;
        elseif shapeCombi == "ball-box" then
          if shape_a == "ball" then
            (gap, contactThreshold, s_rel) := gapBallBox(s_a, geo_a.L, s_b, {geo_b.L, geo_b.W, geo_b.H});
          else
            (gap, contactThreshold, s_rel) := gapBallBox(s_b, geo_b.L, s_a, {geo_a.L, geo_a.W, geo_a.H});
          end if;
        else
          s_rel := s_a - s_b;
          gap := length(s_a - s_b);
          contactThreshold := geo_a.L / 2 + geo_b.L / 2;
        end if;
        L_p := -min(gap - contactThreshold, 0);
        vec_c := normalize(s_rel);
        vec_p := L_p * vec_c;
        contact := if gap - contactThreshold < 0 then true else false;
      end penetration;

      function gapBallBox
      /******************************************************************
      gapBallBox
      ******************************************************************/
        input Real ballCentor[3];
        input Real D "diameter of ball";
        input Real boxCentor[3];
        input Real boxLWH[3];
        output Real gap;
        output Real contactThreshold;
        output SI.Position s_rel[3];
      protected
        Real XYZmin[3];
        Real XYZmax[3];
        Real shortestPointOnBoxWithBall[3];
      algorithm
        XYZmin := boxCentor - boxLWH / 2;
        XYZmax := boxCentor + boxLWH / 2;
        for i in 1:3 loop
          shortestPointOnBoxWithBall[i] := max(XYZmin[i], min(ballCentor[i], XYZmax[i]));
        end for;
        gap := length(ballCentor - shortestPointOnBoxWithBall);
        contactThreshold := D / 2;
        s_rel := shortestPointOnBoxWithBall - ballCentor;
      end gapBallBox;
    equation
/******************************************************************
      following is equation section
      When ball and box contact, ball have to connect contact.falnge_b, and box have to connect contact.flange_a
      ******************************************************************/
      s_rel = s_a - s_b;
      flange_b.f = f;
      flange_a.f = -f;
      (vec_p, contact) = penetration(s_a, s_b, flange_a.geometry, flange_b.geometry);
      v_rel = der(s_rel);
/*
      f_c = smooth(1, noEvent(if contact then c*vec_p else zeros(3)));
      f_d2 = if contact then -d*v_rel else zeros(3);
      f_d = smooth(0, noEvent(if contact then f_d2 else zeros(3)));
      f = f_c + f_d;
      //*/
///* 以下の式は簡略式
      f_c = c * vec_p;
      f_d2 = -d * v_rel;
      f_d = if contact then f_d2 else zeros(3);
      f = f_c + f_d;
//*/
//port
      flange_a.xyz = s_a;
      flange_b.xyz = s_b;
      annotation(
        Icon(graphics = {Polygon(origin = {-1, 9}, fillColor = {255, 255, 0}, fillPattern = FillPattern.Solid, points = {{-55, 87}, {-49, 43}, {-47, 21}, {-95, 15}, {-47, -29}, {-59, -87}, {-13, -37}, {-11, -79}, {5, -47}, {21, -89}, {39, -37}, {87, -65}, {51, -15}, {93, -1}, {53, 15}, {77, 27}, {61, 27}, {85, 79}, {23, 43}, {9, 89}, {-7, 45}, {-55, 87}})}));
    end Contacts;

    model MultiBalls
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
      Interfaces.Flange_a flange_a[n] annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      for i in 1:num loop
        connect(ball[A1[i]].flange_a, contacts[i].flange_a);
        connect(ball[A2[i]].flange_a, contacts[i].flange_b);
      end for;
      for j in 1:n loop
        connect(ball[j].flange_a, flange_a[j]);
      end for;
      annotation(
        Icon(graphics = {Ellipse(origin = {-44, 75}, fillColor = {170, 255, 0}, fillPattern = FillPattern.Sphere, extent = {{-28, 27}, {48, -47}}, endAngle = 360), Ellipse(origin = {76, 17}, fillColor = {0, 255, 255}, fillPattern = FillPattern.Sphere, extent = {{-68, 67}, {32, -31}}, endAngle = 360), Ellipse(origin = {-78, -21}, fillColor = {85, 255, 0}, fillPattern = FillPattern.Sphere, extent = {{-22, 23}, {86, -79}}, endAngle = 360), Ellipse(origin = {36, -25}, fillColor = {255, 0, 0}, fillPattern = FillPattern.Sphere, extent = {{2, -1}, {56, -53}}, endAngle = 360)}, coordinateSystem(initialScale = 0.1)),
        experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002));
    end MultiBalls;

    model FixedMultiBalls
      parameter Integer n = 3;
      parameter Integer num = integer(n * (n - 1) / 2);
      parameter Real position[n, 3] = fill(0.1, n, 3);
      parameter Real v_start[n, 3] = fill(0, n, 3);
      parameter Real D[n] = fill(0.1, n);
      parameter Real m[n] = fill(0.1, n);
      PythagoraDevice.Components.FixedBall ball[n](sphereColor = fill({255, 255, 255}, n), D = D, m = m, position = position, s(fixed = true, start = position), v(fixed = false)) annotation(
        Placement(visible = true, transformation(origin = {-2, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.Flange_a flange_a[n] annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      for j in 1:n loop
        connect(ball[j].flange_a, flange_a[j]);
        ball[j].flange_a.xyz = position[j];
      end for;
      annotation(
        Icon(graphics = {Ellipse(origin = {-44, 75}, fillColor = {85, 0, 0}, fillPattern = FillPattern.Sphere, extent = {{-28, 27}, {48, -47}}, endAngle = 360), Ellipse(origin = {76, 17}, fillPattern = FillPattern.Sphere, extent = {{-68, 67}, {32, -31}}, endAngle = 360), Ellipse(origin = {-78, -21}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Sphere, extent = {{-22, 23}, {86, -79}}, endAngle = 360), Ellipse(origin = {36, -25}, fillColor = {0, 85, 127}, fillPattern = FillPattern.Sphere, extent = {{2, -1}, {56, -53}}, endAngle = 360)}, coordinateSystem(initialScale = 0.1)),
        experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002));
    end FixedMultiBalls;

    model FixedBall
      import Modelica.Mechanics.MultiBody.Types;
      //parameter
      parameter SI.Mass m(min = 0) = 1 "Mass";
      parameter SI.Diameter D(min = 0) = 1 "Diameter";
      parameter Real[3] n = {0, -1, 0} "Direction of movement";
      parameter Real[3] position = {0, 0, 0} "Centor position of object";
      parameter StateSelect stateSelect = StateSelect.default "Priority to use s and v as states" annotation(
        Dialog(tab = "Advanced"));
      parameter Types.Color sphereColor = Modelica.Mechanics.MultiBody.Types.Defaults.BodyColor;
      parameter SI.Force g_f[3] = {0, 0, 0} "Gravity force";
      //physical models
      SI.Velocity v[3](start = {0, 0, 0}) "Absolute velocity of component";
      SI.Acceleration a[3](start = {0, 0, 0}) "Absolute acceleration of component";
      SI.Position s[3](start = {0, 0, 0}) "Absolute position of center of component";
      //Animation
      Modelica.Blocks.Interfaces.RealOutput xyz[3] "Postion of ball";
      PythagoraDevice.Visualizers.Sphere sphere1(sphereColor = sphereColor, sphereDiameter = D) annotation(
        Placement(visible = true, transformation(origin = {14, -54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Interfaces.Flange_a flange_a annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    protected
      constant String selectShape = "ball";
      PythagoraDevice.Interfaces.Geometry geometry(shape = selectShape, L = D, W = D, H = D);
    equation
      v = der(s);
      a = der(v);
      xyz = s + position;
      connect(sphere1.u, xyz);
//port
      flange_a.xyz = xyz;
      flange_a.geometry = geometry;
      annotation(
        Icon(graphics = {Ellipse(fillColor = {170, 255, 255}, fillPattern = FillPattern.Sphere, extent = {{-100, 100}, {100, -100}}, endAngle = 360), Line(origin = {179.403, 0.196862}, points = {{-72, 78}, {-72, -82}}, thickness = 1.5, arrow = {Arrow.None, Arrow.Filled}, arrowSize = 18), Text(origin = {98, -86}, extent = {{-16, 8}, {36, -22}}, textString = "g"), Text(origin = {-2, 117}, extent = {{-68, 17}, {68, -17}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)),
        __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"),
        experiment(StartTime = 0, StopTime = 3, Tolerance = 1e-6, Interval = 0.006));
    end FixedBall;
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

    connector Flange "One-dimensional translational flange"
      SI.Position xyz[3] "Absolute position of center of object";
      flow SI.Force f[3] "Cut force directed into flange";
      Geometry geometry "Geometry infomation of object";
      annotation(
        Documentation(info = "<html>
    <p>
    This is a connector for 1D translational mechanical systems.
    It has no icon definition and is only used by inheritance from
    flange connectors to define different icons.
    </p>
    <p>
    The following variables are defined in this connector:
    </p>
    
    <blockquote><pre>
    s: Absolute position of the flange in [m]. A positive translation
    means that the flange is translated along the flange axis.
    f: Cut-force in direction of the flange axis in [N].
    </pre></blockquote>
    </html>"));
    end Flange;

    connector Flange_a "One-dimensional translational flange (left, flange axis directed INTO cut plane)"
      extends Flange;
      annotation(
        defaultComponentName = "flange_a",
        Icon(coordinateSystem(initialScale = 0.1), graphics = {Rectangle(lineColor = {0, 127, 0}, fillColor = {0, 97, 0}, fillPattern = FillPattern.Solid, extent = {{-100, -100}, {100, 100}})}),
        Diagram(coordinateSystem(initialScale = 0.1), graphics = {Rectangle(lineColor = {0, 127, 0}, fillColor = {0, 127, 0}, fillPattern = FillPattern.Solid, extent = {{-40, -40}, {40, 40}}), Text(lineColor = {0, 127, 0}, extent = {{-160, 110}, {40, 50}}, textString = "%name")}));
    end Flange_a;

    connector Flange_b "One-dimensional translational flange (right, flange axis directed OUT OF cut plane)"
      extends Flange;
      annotation(
        defaultComponentName = "flange_b",
        Icon(coordinateSystem(initialScale = 0.1), graphics = {Rectangle(lineColor = {0, 79, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, -100}, {100, 100}})}),
        Diagram(coordinateSystem(initialScale = 0.1), graphics = {Rectangle(lineColor = {0, 127, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-40, -40}, {40, 40}}), Text(lineColor = {0, 127, 0}, extent = {{-40, 110}, {160, 50}}, textString = "%name")}));
    end Flange_b;

    record Geometry
      String shape "shape of object. ex. ball, box";
      SI.Length L "length";
      SI.Length W "wide";
      SI.Length H "height";
    end Geometry;

    package Partials "Partial models"
      extends Modelica.Icons.BasesPackage;

      partial model TwoPort
        //port
        PythagoraDevice.Interfaces.Flange_a flange_a annotation(
          Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        PythagoraDevice.Interfaces.Flange_b flange_b annotation(
          Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation

      end TwoPort;
    end Partials;
  end Interfaces;

  model Enviroment
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
  end Sources;
  annotation(
    uses(Modelica(version = "3.2.3")),
    Icon(graphics = {Ellipse(fillColor = {165, 244, 121}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}, endAngle = 360), Text(origin = {4, -12}, fillColor = {255, 255, 0}, fillPattern = FillPattern.Solid, extent = {{-83, 71}, {81, -69}}, textString = "ピ", fontName = "UD デジタル 教科書体 NK-B")}, coordinateSystem(initialScale = 0.1)));
end PythagoraDevice;
