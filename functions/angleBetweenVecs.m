
function theta = angleBetweenVecs(v1,v2)
% angleBetweenVecs - Given 2 vectors, finds the angle between them
% theta = angleBetweenVecs(v1,v2)
theta = acos(dot(v1,v2)/(norm(v1)*norm(v2))) * (180/pi);

end

 