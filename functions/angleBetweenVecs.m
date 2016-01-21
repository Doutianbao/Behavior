
function theta = angleBetweenVecs(v1,v2)
theta = acos(dot(v1,v2)/(norm(v1)*norm(v2))) * (180/pi);

end

 