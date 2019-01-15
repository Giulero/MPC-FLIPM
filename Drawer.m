classdef Drawer < handle
  properties
  end
  methods

    function self= Drawer()
    end

    function drawing = drawCircle2D(~, x, y, r,color)
      resolution = 360;
      delta = 2*pi/resolution;

      points = zeros(resolution,2);
      index =1;
      for i = 0.0:delta:2*pi
        p_x= r*cos(i);
        p_y= r*sin(i);
        points(index,:) = [p_x,p_y];
        index = index + 1;
      end
      points(:,1)= points(:,1) + x;
      points(:,2)= points(:,2) + y;
      drawing = plot(points(:,1),points(:,2),'Color',color,'LineWidth', 2);
    end

    function drawing = drawPoligon2D(self, points ,color)
      drawing = [ ];
      for i = 1:size( points,1)-1
        drawing(end)= drawLine2D(self, points(i,:), points(i+1,:), color);
      end
      drawing(end)= drawLine2D(self, points(size( points,1),:), points(1,:), color);
    end

    function drawing = drawRectangle2D(self, points ,color)
      l1= drawLine2D(self, points(1,:), points(2,:), color);
      l2= drawLine2D(self, points(2,:), points(3,:), color);
      l3= drawLine2D(self, points(3,:), points(4,:), color);
      l4= drawLine2D(self, points(4,:), points(1,:), color);
      drawing = [ l1 ; l2 ; l3 ; l4];
    end

    function drawing= drawLine2D(~, first, second, color)
      drawing= line( [ first(1,1), second(1,1)],[ first(1,2), second(1,2)], 'Color', color , 'LineWidth',3);
    end

    function drawing = drawArrow(self,start,edge,headWidth,color)
      headEdgeAngle = 30*pi/180;
      headShaft = headWidth/ ( 2* sin( headEdgeAngle/2));

      lengthArrow = sqrt((edge(1,1)-start(1,1))^2+(edge(1,2)-start(1,2))^2+(edge(1,3)-start(1,3))^2);
      scaleHead = 100*headShaft/lengthArrow;
      scaleTail = 100- scaleHead;
      baseHeadPoint = (start' + scaleTail*( edge'-start')/100)';
      l1 = self.drawLine3D( start,baseHeadPoint,color);
      theta = atan2( edge(1,2)-start(1,2),edge(1,1)-start(1,1));
      sideDisplacement = [cos( theta - pi/2)*headWidth/2; sin(theta - pi/2)*headWidth/2; 0 ];
      rightPoint =  (baseHeadPoint' + sideDisplacement)';
      leftPoint =  (baseHeadPoint' - sideDisplacement)';
      l2 = self.drawLine3D( leftPoint,rightPoint,color);
      l3 = self.drawLine3D( rightPoint,edge,color);
      l4 = self.drawLine3D( edge,leftPoint,color);
      drawing = [ l1 ; l2 ; l3 ; l4];
    end
  end
end

