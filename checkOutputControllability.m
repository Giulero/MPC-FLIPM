function r = checkOutputControllability(A, B, C)
    n_state = size(B,1);
    ctrb = zeros(1, n_state+1);
    for i= 1: (n_state)
        ctrb(:,i) = C*A^(i-1)*B;
    end
    r = rank(ctrb);
end

