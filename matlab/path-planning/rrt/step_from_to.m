function point_ = step_from_to(p1,p2,epsilon_)

if euclidean_dist(p1,p2) < epsilon_
    point_ = p2;
else
    theta_ = atan2(p2(2)-p1(2),p2(1)-p1(1));
    point_(1) = p1(1) + epsilon_*cos(theta_);
    point_(2) = p1(2) + epsilon_*sin(theta_);
end

end