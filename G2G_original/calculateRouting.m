function [adjacencyMat] = calculateRouting(cat, model, surfaceFlow)
    catchmentSize = size(cat.river);
    area = catchmentSize(1) * catchmentSize(2);
    if surfaceFlow
        speed = (1 - cat.river) .* model.cl + cat.river .* model.cr;
    else
        speed = (1 - cat.river) .* model.clb + cat.river .* model.crb;
    end

    flowTime = inf(catchmentSize(1), catchmentSize(2));
    flowTime(cat.outlet(1), cat.outlet(2)) = 0;
    outId = zeros(area, 1);
    
    [~, outId] = examineTile(flowTime, speed, cat.outlet, outId);
    outletId = getId(cat.outlet, catchmentSize);
    otherIds = [1:(outletId-1), (outletId+1):area];
    
    adjacencyMat = sparse(outId(otherIds), otherIds, ones(1, area - 1), ...
                          area, area);

    function [flowTime, outId] = examineTile(flowTime, speed, tile, outId)
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
                outId(getId(neighbour, domainSize)) = getId(tile, domainSize);
                [flowTime, outId] = examineTile(flowTime, speed, ...
                                    [neighbour(1), neighbour(2)], outId);
              end
            end
          end
        end
      end
    end

    function [id] = getId(tile, domainSize)
        id = (tile(2) - 1) * domainSize(1) + tile(1);
    end

end

