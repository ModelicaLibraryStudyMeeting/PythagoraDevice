package PythagoraDeviceTest
  import SI = Modelica.SIunits;
  import Modelica.Math.Vectors.*;
  import Cv = Modelica.SIunits.Conversions;
  import C = Modelica.Constants;
  import Modelica.Mechanics.MultiBody.Frames;
  import Modelica.Mechanics.MultiBody.Types;
  import Modelica.Mechanics.MultiBody.*;
  import Modelica.Constants.*;

  package ComponentTest
    model BallTest
      PythagoraDevice.Components.BallTest1 ball(angles_g(start = {0, 0, eps}), g_f = {0, 0, 0}) annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Sources.InputForceTorque inputForceTorque(Fx = 0, tx = 0, tz = 0) annotation(
        Placement(visible = true, transformation(origin = {-54, -38}, extent = {{-10, -6}, {10, 6}}, rotation = 0)));
      inner Modelica.Mechanics.MultiBody.World world(axisLength = 0.01) annotation(
        Placement(visible = true, transformation(origin = {68, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(inputForceTorque.Center_a, ball.center) annotation(
        Line(points = {{-44, -38}, {0, -38}, {0, 0}}, color = {0, 127, 0}));
      annotation(
        experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.0001));
    end BallTest;

    model BallTest2
      PythagoraDevice.Components.Ball ball(angles_g(displayUnit = "rad", start = {0, 0, 0}), g_f = {0, 0, 0}, w_g_ini = {6.38, -6.81, -25.54}) annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      /*start = {32.7716438460142, -16.0819704180493, 13.8593121283851}*/
      PythagoraDevice.Sources.InputForceTorque inputForceTorque(Fx = 0, tx = 0, tz = 0) annotation(
        Placement(visible = true, transformation(origin = {-54, -38}, extent = {{-10, -6}, {10, 6}}, rotation = 0)));
      inner Modelica.Mechanics.MultiBody.World world(axisLength = 0.01) annotation(
        Placement(visible = true, transformation(origin = {68, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(inputForceTorque.Center_a, ball.center) annotation(
        Line(points = {{-44, -38}, {0, -38}, {0, 0}}, color = {0, 127, 0}));
      annotation(
        experiment(StartTime = 0, StopTime = 0.3, Tolerance = 1e-06, Interval = 0.0001));
    end BallTest2;

    model BoxTest1
      PythagoraDevice.Components.Wall_rotational wall_rotational annotation(
        Placement(visible = true, transformation(origin = {4, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Sources.InputForceTorque inputForceTorque annotation(
        Placement(visible = true, transformation(origin = {-36, -24}, extent = {{-10, -6}, {10, 6}}, rotation = 0)));
      inner Modelica.Mechanics.MultiBody.World world(axisLength = 0.01) annotation(
        Placement(visible = true, transformation(origin = {68, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(inputForceTorque.Center_a, wall_rotational.Center_a) annotation(
        Line(points = {{-26, -24}, {4, -24}, {4, 2}}, color = {0, 127, 0}));
      annotation(
        experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002));
    end BoxTest1;

    package RotationTest
      model RotBox
        parameter Real r_shape[3] = {-length / 2, 0, 0} "set the centor to origin";
        constant String selectShape = "box";
        parameter SI.Length length = 0.05;
        parameter SI.Length width = 0.01;
        parameter SI.Length height = 0.025;
        //  parameter Real angles_g_start[3]={32.7716438460142, -16.0819704180493, 13.8593121283851};
        parameter Real angles_g_start[3] = {0, 3.141592 / 2, 0};
        parameter Real changeTime = 1;
        parameter SI.Position position_ini[3] = {0, 0, 0} "centor of box";
        parameter Real color[3] = {255, 255, 0} annotation(
          Dialog(tab = "Visualization"));
        Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape vis(shapeType = selectShape, color = color, length = length, width = width, height = height, lengthDirection = {1, 0, 0}, widthDirection = {0, 1, 0}, r_shape = r_shape, r = position_ini, R = R);
        Frames.Orientation R;
        SI.Angle angles_g[3](start = angles_g_start);
        SI.AngularVelocity w_g[3] "Absolute angular velocity of frame_a resolved in world fram";
        Frames.Orientation R21;
        Frames.Orientation R_a;
        Real angles_a;
      equation
//  w_g = if time < changeTime then {0, 0, time} else {0, time, 0};
        w_g = {0, 0, time};
        w_g = der(angles_g);
        R_a = Frames.axesRotations(sequence = {1, 2, 3}, angles = -angles_g, der_angles = w_g);
        angles_a = Frames.resolve1(R_a, angles_g);
        R_b = Frames.axesRotations(sequence = {3, 2, 1}, angles = angles_a, der_angles = w_g);
        w1 = {0, 0, time};
        w1 = der(angles1);
        angles0_1 = angles0_0 + resolve(R, angles1);
        R = Frames.axesRotations(sequence = {1, 2, 3}, angles = angles_g0_0, der_angles = w_g);
        R2 = Frames.axesRotations(sequence = {1, 2, 3}, angles = angles_g0_1, der_angles = w_g);
        annotation(
          experiment(StartTime = 0, StopTime = 2, Tolerance = 1e-06, Interval = 0.004));
      end RotBox;

      model RotBox2
        parameter Real r_shape[3] = {-length / 2, 0, 0} "set the centor to origin";
        constant String selectShape = "box";
        parameter SI.Length length = 0.05;
        parameter SI.Length width = 0.01;
        parameter SI.Length height = 0.025;
        //  parameter Real angles_g_start[3]={32.7716438460142, -16.0819704180493, 13.8593121283851};
        parameter Real angles_g_start[3] = {0, 3.141592 / 2, 0};
        parameter Real changeTime = 1;
        parameter SI.Position position_ini[3] = {0, 0, 0} "centor of box";
        parameter Real color[3] = {255, 255, 0} annotation(
          Dialog(tab = "Visualization"));
        Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape vis(shapeType = selectShape, color = color, length = length, width = width, height = height, lengthDirection = {1, 0, 0}, widthDirection = {0, 1, 0}, r_shape = r_shape, r = position_ini, R = R2);
        Frames.Orientation R1;
        Frames.Orientation R2;
        SI.Angle angles1[3];
        SI.Angle angles2[3];
        SI.Angle angles_g0_1[3];
        parameter SI.Angle angles_g0_0[3] = {0, 3.141592 / 2, 0};
        SI.AngularVelocity w1[3] "Absolute angular velocity of frame_a resolved in world fram";
      equation
        w1 = {0, 0, time};
        w1 = der(angles1);
        R1 = Frames.axesRotations(sequence = {1, 2, 3}, angles = angles_g0_0, der_angles = w1);
        angles2 = Frames.resolve1(R1, angles1);
        angles_g0_1 = (+angles_g0_0) + Frames.resolve1(R1, angles1);
        R2 = Frames.axesRotations(sequence = {1, 2, 3}, angles = angles_g0_1, der_angles = w1);
        annotation(
          experiment(StartTime = 0, StopTime = 2, Tolerance = 1e-06, Interval = 0.004));
      end RotBox2;

      model RotBox3
        import Modelica.Math.Vectors.*;
        parameter Real r_shape[3] = {-length / 2, 0, 0} "set the centor to origin";
        constant String selectShape = "box";
        parameter SI.Length length = 0.05;
        parameter SI.Length width = 0.01;
        parameter SI.Length height = 0.025;
        //  parameter Real angles_g_start[3]={32.7716438460142, -16.0819704180493, 13.8593121283851};
        parameter Real angles_g_start[3] = {0, 3.141592 / 2, 0};
        parameter Real changeTime = 1;
        parameter SI.Position position_ini[3] = {0, 0, 0} "centor of box";
        parameter Real color[3] = {255, 255, 0} annotation(
          Dialog(tab = "Visualization"));
        Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape vis(shapeType = selectShape, color = color, length = length, width = width, height = height, lengthDirection = {1, 0, 0}, widthDirection = {0, 1, 0}, r_shape = r_shape, r = position_ini, R = R);
        Frames.Orientation R;
        SI.Angle angles1[3](start = {0, 0, 0});
        SI.Angle angles_g[3];
        parameter SI.Angle angles_g0[3] = {0, 3.141592 / 2, 0};
        SI.AngularVelocity w1[3] "Absolute angular velocity of frame_a resolved in world fram";
        Real Rx[3, 3];
        Real Ry[3, 3];
        Real Ry1[3, 3];
        Real Rz[3, 3];
        Modelica.Mechanics.MultiBody.Frames.Quaternions.Orientation Q;
      equation
        Ry1 = [Modelica.Math.cos(-angles_g0[2]), 0, Modelica.Math.sin(-angles_g0[2]); 0, 1, 0; -Modelica.Math.sin(-angles_g0[2]), 0, Modelica.Math.cos(-angles_g0[2])];
        Rx = [1, 0, 0; 0, Modelica.Math.cos(angles_g[1]), -Modelica.Math.sin(angles_g[1]); 0, Modelica.Math.sin(angles_g[1]), Modelica.Math.cos(angles_g[1])];
        Ry = [Modelica.Math.cos(angles_g[2]), 0, Modelica.Math.sin(angles_g[2]); 0, 1, 0; -Modelica.Math.sin(angles_g[2]), 0, Modelica.Math.cos(angles_g[2])];
        Rz = [Modelica.Math.cos(angles_g[3]), -Modelica.Math.sin(angles_g[3]), 0; Modelica.Math.sin(angles_g[3]), Modelica.Math.cos(angles_g[3]), 0; 0, 0, 1];
        w1 = {0, 0, 3.14 / 10};
        w1 = der(angles1);
        angles_g = Ry * angles1;
// R = Frames.axesRotations(sequence = {3, 2, 1}, angles = angles_g0 + angles1, der_angles = w1);
        Q = Modelica.Mechanics.MultiBody.Frames.Quaternions.planarRotation(normalize(angles_g0 + angles1), Modelica.Math.Vectors.length(angles_g0 + angles1));
        R = Frames.from_Q(Q, w1);
//R.T=Rz*Ry*Rx; //固定角 x-y-zの順で回転
//    R.T=Rx*Ry*Rz; //オイラー角 x-y-zの順で回転
//    R.w=w1;
        annotation(
          experiment(StartTime = 0, StopTime = 2, Tolerance = 1e-06, Interval = 0.004));
      end RotBox3;

      model RotBox4
        import Modelica.Math.Vectors.*;
        parameter Real r_shape[3] = {-length / 2, 0, 0} "set the centor to origin";
        constant String selectShape = "box";
        parameter SI.Length length = 0.05;
        parameter SI.Length width = 0.01;
        parameter SI.Length height = 0.025;
        //  parameter Real angles_g_start[3]={32.7716438460142, -16.0819704180493, 13.8593121283851};
        parameter Real angles_g_start[3] = {0, 3.141592 / 2, 0};
        parameter Real changeTime = 1;
        parameter SI.Position position_ini[3] = {0, 0, 0} "centor of box";
        parameter Real color[3] = {255, 255, 0} annotation(
          Dialog(tab = "Visualization"));
        Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape vis(shapeType = selectShape, color = color, length = length, width = width, height = height, lengthDirection = {1, 0, 0}, widthDirection = {0, 1, 0}, r_shape = r_shape, r = position_ini, R = R);
        Frames.Orientation R;
        SI.Angle angles_g[3](start = {0, 0, 0});
        SI.Angle d_angles_g[3](start = {0, 0, 0});
        //  Real delay_angles_g[3];
        parameter Real dt = 0.01;
        parameter SI.Angle angles_g0[3] = {0, 3.141592 / 2, 0};
        SI.AngularVelocity w1[3] "Absolute angular velocity of frame_a resolved in world fram";
        Real Rx[3, 3];
        Real Ry[3, 3];
        Real Ry1[3, 3];
        Real Rz[3, 3];
        //  Modelica.Mechanics.MultiBody.Frames.Quaternions.Orientation Q;
      equation
        Ry1 = [Modelica.Math.cos(angles_g0[2]), 0, Modelica.Math.sin(angles_g0[2]); 0, 1, 0; Modelica.Math.sin(angles_g0[2]), 0, Modelica.Math.cos(angles_g0[2])];
        Rx = [1, 0, 0; 0, Modelica.Math.cos(d_angles_g[1]), Modelica.Math.sin(d_angles_g[1]); 0, -Modelica.Math.sin(d_angles_g[1]), Modelica.Math.cos(d_angles_g[1])];
        Ry = [Modelica.Math.cos(d_angles_g[2]), 0, -Modelica.Math.sin(d_angles_g[2]); 0, 1, 0; Modelica.Math.sin(d_angles_g[2]), 0, Modelica.Math.cos(d_angles_g[2])];
        Rz = [Modelica.Math.cos(d_angles_g[3]), Modelica.Math.sin(d_angles_g[3]), 0; -Modelica.Math.sin(d_angles_g[3]), Modelica.Math.cos(d_angles_g[3]), 0; 0, 0, 1];
        w1 = if time < 1 then {0, 3.14 / 2, 0} else {0, 0, 3.14 / 2};
        w1 = der(angles_g);
        for i in 1:3 loop
          d_angles_g[i] = delay(w1[i], dt) * dt;
        end for;
        R.T = Ry1 * Rz;
/* 
           固定座標系y軸に90°回転した後に、固定座標系z軸にゆっくり回転　回転させるのを右からかける？
           じゃ、初期回転angles0はR0.T、後の回転angles1はR1.Tとすると座標変換行列R.Tは
           R.T = R0.T * R1.T
           問題はR0.Tはどの順番で回せばいいのか？
           どの順番でかければいいのか？
        */
//    R.T = Rz * Ry * Rx; 固定角 x-y-zの順で回転
//    R.T=Rx*Ry*Rz; //オイラー角 x-y-zの順で回転
        R.w = w1;
        annotation(
          experiment(StartTime = 0, StopTime = 2, Tolerance = 1e-06, Interval = 0.002));
      end RotBox4;

      model RotBox5
        import Modelica.Math.Vectors.*;
        parameter Real r_shape[3] = {-length / 2, 0, 0} "set the centor to origin";
        constant String selectShape = "box";
        parameter SI.Length length = 0.05;
        parameter SI.Length width = 0.01;
        parameter SI.Length height = 0.025;
        //  parameter Real angles_g_start[3]={32.7716438460142, -16.0819704180493, 13.8593121283851};
        parameter Real angles_g_start[3] = {0, 3.141592 / 2, 0};
        parameter Real changeTime = 1;
        parameter SI.Position position_ini[3] = {0, 0, 0} "centor of box";
        parameter Real color[3] = {255, 255, 0} annotation(
          Dialog(tab = "Visualization"));
        Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape vis(shapeType = selectShape, color = color, length = length, width = width, height = height, lengthDirection = {1, 0, 0}, widthDirection = {0, 1, 0}, r_shape = r_shape, r = position_ini, R = R);
        Frames.Orientation R;
        Frames.Orientation d_R;
        //  SI.Angle angles_g[3](start = {0, 0, 0});
        SI.Angle d_angles_g[3](start = {0, 0, 0});
        //  Real delay_angles_g[3];
        parameter Real dt = 0.01;
        parameter SI.Angle angles_g0[3] = {0, 3.141592 / 2, 0};
        SI.AngularVelocity w1[3] "Absolute angular velocity of frame_a resolved in world fram";
        Real Rx[3, 3];
        Real Ry[3, 3];
        Real Ry1[3, 3];
        Real Rz[3, 3];
        Real T_delay[3, 3];
        //  Modelica.Mechanics.MultiBody.Frames.Quaternions.Orientation Q;
      initial equation
        R = Frames.axesRotations({1, 2, 3}, zeros(3), zeros(3));
      equation
        Ry1 = [Modelica.Math.cos(angles_g0[2]), 0, Modelica.Math.sin(angles_g0[2]); 0, 1, 0; Modelica.Math.sin(angles_g0[2]), 0, Modelica.Math.cos(angles_g0[2])];
        Rx = [1, 0, 0; 0, Modelica.Math.cos(d_angles_g[1]), Modelica.Math.sin(d_angles_g[1]); 0, -Modelica.Math.sin(d_angles_g[1]), Modelica.Math.cos(d_angles_g[1])];
        Ry = [Modelica.Math.cos(d_angles_g[2]), 0, -Modelica.Math.sin(d_angles_g[2]); 0, 1, 0; Modelica.Math.sin(d_angles_g[2]), 0, Modelica.Math.cos(d_angles_g[2])];
        Rz = [Modelica.Math.cos(d_angles_g[3]), Modelica.Math.sin(d_angles_g[3]), 0; -Modelica.Math.sin(d_angles_g[3]), Modelica.Math.cos(d_angles_g[3]), 0; 0, 0, 1];
        w1 = if time < 1 then {0, 3.14 / 2, 0} else {0, 0, 3.14 / 2};
//  w1 = der(angles_g);
        for i in 1:3 loop
          d_angles_g[i] = delay(w1[i], dt) * dt;
          for j in 1:3 loop
            T_delay[i, j] = delay(R.T[i, j], dt);
          end for;
        end for;
        d_R = Frames.axesRotations({1, 2, 3}, d_angles_g, w1);
        R.T = T_delay * d_R.T;
        R.w = w1;
//  R.T = Ry1 * Rz;
/* 
           固定座標系y軸に90°回転した後に、固定座標系z軸にゆっくり回転　回転させるのを右からかける？
           じゃ、初期回転angles0はR0.T、後の回転angles1はR1.Tとすると座標変換行列R.Tは
           R.T = R0.T * R1.T
           問題はR0.Tはどの順番で回せばいいのか？
           どの順番でかければいいのか？
        */
//    R.T = Rz * Ry * Rx; 固定角 x-y-zの順で回転
//    R.T=Rx*Ry*Rz; //オイラー角 x-y-zの順で回転
        annotation(
          experiment(StartTime = 0, StopTime = 2, Tolerance = 1e-06, Interval = 0.01));
      end RotBox5;

      model MSLBoxTest1
        inner Modelica.Mechanics.MultiBody.World world annotation(
          Placement(visible = true, transformation(origin = {-76, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Mechanics.MultiBody.Parts.BodyBox bodyBox(r = {1, 0, 0}) annotation(
          Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation

        annotation(
          experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002));
      end MSLBoxTest1;

      model GetPreStepValue
      
      Real a;
      Real a_pre(start=1);
      parameter Real b=2;
      
      algorithm
        when sample(0, 0.5) then
          a:=a_pre*b;
          a_pre:=a;
        end when;
      end GetPreStepValue;
      
      model GetPreStepValue2
      
      Real a;
      Real a_pre(start=1);
      Real db(start=0.25);
      Real b_pre;
      Real b;
      Real pre_b;
      Real diff_b;
      
        
      equation
        b=time;
//  pre_b=noEvent(pre(b));
//  diff_b=noEvent(b-pre_b);
      algorithm
        
        when sample(0, 0.25) then
          db := b - b_pre;
          a := a_pre*db;
          a_pre:=a;
          b_pre:=b;
        end when;
        
        pre_b:=noEvent(pre(b));
        diff_b:=noEvent(b-pre_b);
        
      end GetPreStepValue2;
      
      model GetPreStepValue3
        Clock clk1 = Clock(0.1);
        Real a;
        //Real a_pre(start=1);
        //Real diff_a;
        Real y;
        parameter Real b=2;
      
      
      equation
        y=time;
      
        when Clock(0.1) then
          a = previous(y);
        end when;
          
        
      end GetPreStepValue3;
      
      model GetPreStepValue4
      //  Clock clk1 = Clock(0.1);
        Real a;
        //Real a_pre(start=1);
        //Real diff_a;
        Real y, yy;
        Real b;
      
      
      equation
        b=time;
        der(y)=b;
//  when Clock(0.1) then
//    a=0.1*b;
//  end when;
        when sample(0, 0.1) then
          a=0.1*b;
        end when; 
      
      algorithm
        
      
      end GetPreStepValue4;

      model NoClockVsSampleHold
      Clock clk1 = Clock(0.1);
        Clock clk2 = subSample(clk1, 2);
        Real x(start = 0), y(start = 0), z(start = 0);
      equation
        when clk1 then
          x = previous(x) + 0.1;
        end when;
        when clk2 then
          y = noClock(x);
// most recent value of x
          z = sample(hold(x));
// left limit of x (infinitesimally delayed)!
        end when;
      end NoClockVsSampleHold;
      
      model RotBox6
        import Modelica.Math.Vectors.*;
        parameter Real r_shape[3] = {-length / 2, 0, 0} "set the centor to origin";
        constant String selectShape = "box";
        parameter SI.Length length = 0.05;
        parameter SI.Length width = 0.01;
        parameter SI.Length height = 0.025;
      
      
        parameter SI.Position position_ini[3] = {0, 0, 0} "centor of box";
        parameter Real color[3] = {255, 255, 0} annotation(
          Dialog(tab = "Visualization"));
        Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape vis(shapeType = selectShape, color = color, length = length, width = width, height = height, lengthDirection = {1, 0, 0}, widthDirection = {0, 1, 0}, r_shape = r_shape, r = position_ini, R = R);
        Frames.Orientation R;
      //  Real delay_angles_g[3];
        parameter Real dt = 0.01;
        SI.Angle angles1[3](start={0, 0, 0});
        SI.Angle angles2[3](start={0, 0, 0});
        SI.AngularVelocity w1[3] "Absolute angular velocity of frame_a resolved in world fram";
        SI.AngularVelocity w2[3] "Absolute angular velocity of frame_a resolved in world fram";
        Real Ry1[3, 3];
        Real Rz2[3, 3];
        //  Modelica.Mechanics.MultiBody.Frames.Quaternions.Orientation Q;
      equation
//    R.T = Rz * Ry * Rx; 固定角 x-y-zの順で回転
//    R.T=Rx*Ry*Rz; //オイラー角 x-y-zの順で回転
        w1 = {0, 3.14 / 2, 0};
        w2={0, 0, 3.14/2};
        Ry1 = [Modelica.Math.cos(angles1[2]), 0, -Modelica.Math.sin(angles1[2]); 0, 1, 0; Modelica.Math.sin(angles1[2]), 0, Modelica.Math.cos(angles1[2])];
        Rz2 = [Modelica.Math.cos(angles2[3]), Modelica.Math.sin(angles2[3]), 0; -Modelica.Math.sin(angles2[3]), Modelica.Math.cos(angles2[3]), 0; 0, 0, 1];
  if time < 1 then
          w1 = der(angles1);
          R.T = Ry1;
//dummy
          angles2 = zeros(3);
          R.w = w1;
        else
          w2 = der(angles2);
          R.T = Ry1 * Rz2;
//dummy
          angles1 = w1;
          R.w = w2;
        end if;
        annotation(
          experiment(StartTime = 0, StopTime = 2, Tolerance = 1e-06, Interval = 0.002));
      end RotBox6;
      
      model RotBox7
        import Modelica.Math.Vectors.*;
        parameter Real r_shape[3] = {-length / 2, 0, 0} "set the centor to origin";
        constant String selectShape = "box";
        parameter SI.Length length = 0.05;
        parameter SI.Length width = 0.01;
        parameter SI.Length height = 0.025;
      
      
        parameter SI.Position position_ini[3] = {0, 0, 0} "centor of box";
        parameter Real color[3] = {255, 255, 0} annotation(
          Dialog(tab = "Visualization"));
        Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape vis(shapeType = selectShape, color = color, length = length, width = width, height = height, lengthDirection = {1, 0, 0}, widthDirection = {0, 1, 0}, r_shape = r_shape, r = position_ini, R = R);
        Frames.Orientation R;
      
        Frames.Orientation Rpre;
        Frames.Orientation Rcur;
        parameter Real dt = 0.01;
        SI.Angle d_angles[3];
        SI.AngularVelocity wcur[3] "Absolute angular velocity of frame_a resolved in world fram";
      
      initial equation
//  R=Frames.axesRotations(sequence = {1, 2, 3}, angles = {1e-6,1e-6,1e-6}, der_angles = zeros(3));
      equation
      
        R.T=Rpre.T*Rcur.T;
        R.w=wcur;
      
        if time < dt then
          Rpre=Frames.axesRotations(sequence = {1, 2, 3}, angles = {1e-6,1e-6,1e-6}, der_angles = zeros(3));
        else
          for i in 1:3 loop
            for j in 1:3 loop
              Rpre.T[i, j] = delay(R.T[i, j], dt);
            end for;
          end for;
          Rpre.w=wcur;
        end if;
      
        Rcur=Frames.axesRotations(sequence = {1, 2, 3}, angles = d_angles, der_angles = wcur);
        d_angles = wcur*dt;
        if time < 1 then
          wcur={0, 3.14/2, 0};
        else
          wcur={0, 0, 3.14/2};
        end if;
        
      
        annotation(
          experiment(StartTime = 0, StopTime = 2, Tolerance = 1e-06, Interval = 0.01));
      end RotBox7;
      
      model RotBox8
        import Modelica.Math.Vectors.*;
        parameter Real r_shape[3] = {-length / 2, 0, 0} "set the centor to origin";
        constant String selectShape = "box";
        parameter SI.Length length = 0.05;
        parameter SI.Length width = 0.01;
        parameter SI.Length height = 0.025;
      
      
        parameter SI.Position position_ini[3] = {0, 0, 0} "centor of box";
        parameter Real color[3] = {255, 255, 0} annotation(
          Dialog(tab = "Visualization"));
        Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape vis(shapeType = selectShape, color = color, length = length, width = width, height = height, lengthDirection = {1, 0, 0}, widthDirection = {0, 1, 0}, r_shape = r_shape, r = position_ini, R = R);
        Frames.Orientation R;
      
        Frames.Orientation Rpre;
        Frames.Orientation Rcur;
        parameter Real dt = 0.01;
        SI.Angle d_angles[3];
        SI.AngularVelocity wcur[3] "Absolute angular velocity of frame_a resolved in world fram";
      
      initial equation
//  R=Frames.axesRotations(sequence = {1, 2, 3}, angles = {1e-6,1e-6,1e-6}, der_angles = zeros(3));
      equation
      
      
        R.w=wcur;
        
        d_angles = wcur*dt;
        if time < 1 then
          wcur={0, 3.14/2, 0};
        else
          wcur={0, 0, 3.14/2};
        end if;
      
        Rcur = Frames.axesRotations(sequence = {1, 2, 3}, angles = d_angles, der_angles = wcur);
        
      algorithm
        when sample(0, dt) then
          
          if time < dt then
            Rpre := Frames.axesRotations(sequence = {1, 2, 3}, angles = {1e-6, 1e-6, 1e-6}, der_angles = wcur);
          end if;
          R.T := Rpre.T * Rcur.T;
          Rpre.T := R.T;
          Rpre.w := wcur;
          
        end when;  
      
        annotation(
          experiment(StartTime = 0, StopTime = 2, Tolerance = 1e-06, Interval = 0.01));
      end RotBox8;
      
      model RotBox9 "copy from RotBox7"
        import Modelica.Math.Vectors.*;
        parameter Real r_shape[3] = {-length / 2, 0, 0} "set the centor to origin";
        constant String selectShape = "box";
        parameter SI.Length length = 0.05;
        parameter SI.Length width = 0.01;
        parameter SI.Length height = 0.025;
      
      
        parameter SI.Position position_ini[3] = {0, 0, 0} "centor of box";
        parameter Real color[3] = {255, 255, 0} annotation(
          Dialog(tab = "Visualization"));
        Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape vis(shapeType = selectShape, color = color, length = length, width = width, height = height, lengthDirection = {1, 0, 0}, widthDirection = {0, 1, 0}, r_shape = r_shape, r = position_ini, R = R);
        Frames.Orientation R;
      
        Frames.Orientation Rpre;
        Frames.Orientation Rcur;
        parameter Real dt = 0.01;
        SI.Angle d_angles[3];
        SI.AngularVelocity wcur[3] "Absolute angular velocity of frame_a resolved in world fram";
      
        Modelica.Mechanics.MultiBody.Frames.Quaternions.Orientation Qpre;
      //  Modelica.Mechanics.MultiBody.Frames.Quaternions.Orientation Qcur;
      
      equation
      
        R.T=Rpre.T*Rcur.T;
        R.w=wcur;
      
      
     Qpre = Modelica.Mechanics.MultiBody.Frames.Quaternions.planarRotation(normalize({1e-6,1e-6,1e-6}), Modelica.Math.Vectors.length({1e-6,1e-6,1e-6}));
        if time < dt then
      
          Rpre = Frames.from_Q(Qpre, wcur);
      //    Rpre=Frames.axesRotations(sequence = {1, 2, 3}, angles = {1e-6,1e-6,1e-6}, der_angles = zeros(3));
        else
      
          for i in 1:3 loop
            for j in 1:3 loop
              Rpre.T[i, j] = delay(R.T[i, j], dt);
            end for;
          end for;
          Rpre.w=wcur;
        end if;
      
        Rcur=Frames.axesRotations(sequence = {1, 2, 3}, angles = d_angles, der_angles = wcur);
        d_angles = wcur*dt;
        if time < 1 then
          wcur={0, 3.14/2, 0};
        else
          wcur={0, 0, 3.14/2};
        end if;
        
      
        annotation(
          experiment(StartTime = 0, StopTime = 2, Tolerance = 1e-06, Interval = 0.01));
      end RotBox9;
    end RotationTest;
  end ComponentTest;

  package ContactTest
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

    model FallingBallOntoBox2 "Falling the ball onto a plane"
      //  inner Modelica.Mechanics.MultiBody.World world annotation(
      //    Placement(visible = true, transformation(origin = {-80, 78}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.FixedBox fixedBoxB(height = 0.05, length = 0.2, position_ini = {0, -0.005, 0}, width = 0.01) annotation(
        Placement(visible = true, transformation(origin = {24, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.ContactsBallFixedBox contactsBallFixedBoxB annotation(
        Placement(visible = true, transformation(origin = {98, -14}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      PythagoraDevice.Components.Ball ball_B(D = 11.3e-3, eps = 7, m = 1.89e-3, position_ini = {0, 0.1, 0}, v_g(fixed = true), v_g_ini = {0, 0, 0}, w_g_ini = {0, 0, -30}) annotation(
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
        experiment(StartTime = 0, StopTime = 0.15, Tolerance = 1e-06, Interval = 1e-05),
        Diagram(graphics = {Rectangle(extent = {{28, -12}, {28, -12}}), Text(origin = {-44, 49}, extent = {{-44, 5}, {44, -5}}, textString = "With rotation"), Text(origin = {98, 49}, extent = {{-44, 5}, {44, -5}}, textString = "No rotation")}, coordinateSystem(extent = {{-200, -100}, {200, 100}})),
        Icon(coordinateSystem(extent = {{-200, -100}, {200, 100}})));
    end FallingBallOntoBox2;

    model RollingBallOnBox "Falling the ball onto a plane"
      //  inner Modelica.Mechanics.MultiBody.World world annotation(
      //    Placement(visible = true, transformation(origin = {-80, 78}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.FixedBox fixedBoxB(height = 0.3, length = 0.3, position_ini = {0, -0.011, 0}, width = 0.01) annotation(
        Placement(visible = true, transformation(origin = {-44, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.Ball ball_A(D = 11.3e-3, color1 = {0, 255, 255}, m = 1.89e-3, position_ini = {0, 0, 0}, v_g(fixed = true), v_g_ini = {0.10, 0, 0.10}, w_g_ini = {0, 0, 0}) annotation(
        Placement(visible = true, transformation(origin = {-44, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PythagoraDevice.Components.ContactsBallFixedBox contactsBallFixedBoxA annotation(
        Placement(visible = true, transformation(origin = {-44, -14}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    equation
      connect(ball_A.center, contactsBallFixedBoxA.center_a) annotation(
        Line(points = {{-44, 28}, {-44, -4}}, color = {0, 127, 0}));
      connect(contactsBallFixedBoxA.center_b, fixedBoxB.center) annotation(
        Line(points = {{-44, -24}, {-44, -50}}, color = {0, 79, 0}));
    protected
      annotation(
        experiment(StartTime = 0, StopTime = 0.16, Tolerance = 1e-06, Interval = 0.0001),
        Diagram(graphics = {Rectangle(extent = {{28, -12}, {28, -12}}), Text(origin = {-44, 49}, extent = {{-44, 5}, {44, -5}}, textString = "With rotation")}, coordinateSystem(extent = {{-200, -100}, {200, 100}})),
        Icon(coordinateSystem(extent = {{-200, -100}, {200, 100}})));
    end RollingBallOnBox;

  end ContactTest;
  annotation(
    uses(Modelica(version = "3.2.3")));
end PythagoraDeviceTest;
