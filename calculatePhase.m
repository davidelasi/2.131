function phase = calculatePhase(dx,dy)
   phase = zeros(length(dx),1);
   for i = 2:length(dx)
       phase(i) = atan2(dy(i), dx(i));
   end
end