function ms0 = irc_mirror_mode_replace_atoms(ms0,replacelist)

    for aind=1:numel(replacelist.A)
        buf = [ ms0.x(replacelist.A(aind)) ms0.y(replacelist.A(aind)) ms0.z(replacelist.A(aind))];
        ms0.x(replacelist.A(aind)) = ms0.x(replacelist.B(aind));
        ms0.y(replacelist.A(aind)) = ms0.y(replacelist.B(aind));
        ms0.z(replacelist.A(aind)) = ms0.z(replacelist.B(aind));
        ms0.x(replacelist.B(aind)) = buf(1);
        ms0.y(replacelist.B(aind)) = buf(2);
        ms0.z(replacelist.B(aind)) = buf(3);
    end

end
