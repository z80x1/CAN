%garb
function rs = check_rect_overlap(rect1, rect2)

%checks if two rectangles overlap
%rect = [xstart ystart xsize ysize]

for i=1:numel()
                e2_x1 = e2(1); e2_x2 = e2(1)+e2(3);
                e2_y1 = e2(2); e2_y2 = e2(2)+e2(4);
                if e2_x1 > exts(it).xstart && e2_x1 < exts(it).xstart
