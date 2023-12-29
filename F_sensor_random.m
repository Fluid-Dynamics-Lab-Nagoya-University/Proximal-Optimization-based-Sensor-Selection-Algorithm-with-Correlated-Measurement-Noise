function [isensors]=F_sensor_random( U, p, s)
[n,~]=size(U);
isensors=[];

for i=1:p
   iflag=1;
   while iflag==1
      iflag=0;
      ii=int16(rand*n/s);      
      for iii=1:i-1
         if isensors(iii)==ii
            iflag=1;
         end
      end
      if ii==0
          iflag=1;
      end
   end
   isensors=[isensors ; ii];
   
end