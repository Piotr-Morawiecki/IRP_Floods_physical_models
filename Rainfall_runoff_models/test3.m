domainSize = [20, 20];
conductivity = ones(domainSize(1), domainSize(2));
conductivity(3:16, 4:8) = 0.2;
conductivity(3:5, 4:13) = 0.2;
speed = 1./conductivity;
direction = zeros(domainSize(1), domainSize(2), 2);
flowTime = inf(domainSize(1), domainSize(2));

sink = [20, 15];
flowTime(sink(1), sink(2)) = 0;

[flowTime, direction] = examineTile(flowTime, direction, speed, sink);

subplot(2,2,1)
heatmap(1:20, flip(1:20), flip(conductivity,1));
subplot(2,2,2)
heatmap(1:20, flip(1:20), flip(flowTime,1));
subplot(2,2,3)
quiver(direction(:,:,2),direction(:,:,1));
axis([0 21 0 21])

function [flowTime, flowDir] = examineTile(flowTime, flowDir, speed, tile)
  domainSize = size(flowTime);
  for i = [-1 0 1]
    for j = [-1 0 1]
      if (i~=0 || j~=0)
        distance = sqrt(i^2 + j^2);
        neighbour = tile + [i j];
        if (neighbour(1) > 0 && neighbour(1) <= domainSize(1) && ...
            neighbour(2) > 0 && neighbour(2) <= domainSize(2))
          newTime = distance * ( speed(neighbour(1), neighbour(2))/2 + ...
              speed(tile(1), tile(2))/2 ) + flowTime(tile(1), tile(2));
          if ( newTime < flowTime(neighbour(1), neighbour(2)))
            flowTime(neighbour(1), neighbour(2)) = newTime;
            flowDir(neighbour(1), neighbour(2), :) = [-i, -j]; 
            [flowTime, flowDir] = examineTile(flowTime, flowDir, speed, [neighbour(1), neighbour(2)]);
          end
        end
      end
    end
  end
end
