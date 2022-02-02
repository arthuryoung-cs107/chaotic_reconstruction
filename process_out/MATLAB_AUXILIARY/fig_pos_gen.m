function fig_pos_out = fig_pos_gen(rows_, cols_)
  rows = rows_;
  if rows > 10;
    rows = 10;
  end
  cols = cols_;
  if cols > 10;
    cols = 10;
  end

  fig_pos_out = zeros(rows*cols, 4);

  H = 850;
  W = 1500;
  del = 10;
  hbar = 80;

  h = floor(H/rows-hbar);
  w = floor((W)/cols);

  for i=1:rows*cols
    fig_pos_out(i, 3) = w;
    fig_pos_out(i, 4) = h;
  end
  k = 1;
  p = rows*cols;
  for i=1:rows
    for j=1:cols
      fig_pos_out(k, 1) = w*(j-1);
      fig_pos_out(p, 2) = (h + hbar)*(i-1);
      k = k + 1;
      p = p - 1;
    end
  end

end
