% 生成一条过p1和p2的直线Ax+By+C=0，并且方向是从y值小的点指向y值大的点，返回系数[A; B; C]
function ABC = line_by_2p(p1, p2)
    ABC = [p2(2)-p1(2); p1(1)-p2(1); p2(1)*p1(2)-p1(1)*p2(2)];
    if p1(2) > p2(2)
        ABC = -ABC;
    end
end