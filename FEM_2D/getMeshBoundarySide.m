function S = getMeshBoundarySide(M, i)

%Index + Physical mark
S.idxA = M.LINES(i, 1);
S.idxB = M.LINES(i, 2);
S.mark = M.LINES(i, 3);

%Pocatecni + koncovy bod usecky
S.A = M.POS(S.idxA, 1:2);
S.B = M.POS(S.idxB, 1:2);

%Stred usecky
S.MidP = (S.A + S.B) / 2;

%Length usecky
S.ds = hypot(S.B(1) - S.A(1), S.B(2) - S.A(2));

% Normala - jednotkovy vekor vnejsi normaly ==>> {normala / vzdalenost} = jednotkovy vektor;
S.nn = [ S.B(2) - S.A(2), S.A(1) - S.B(1) ] / S.ds;

end