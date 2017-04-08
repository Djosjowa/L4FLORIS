
function globalOffset = write_compact_lifetime_damage_table( lifetimeType, globalOffset, sheetObj, FileInfo, Fatigue, realFmt )

   nFiles    = size( FileInfo.FileName, 1 );
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Create table headers, etc.
  
  
   for iGroup=1:Fatigue.nGroups
      nGroupChannels = length(Fatigue.Groups(iGroup).channelIndices);
      
      switch lifetimeType
         case 1
            write_excel_cells( sheetObj, sprintf('%s Lifetime Damage (-) for various S/N Curves', char(Fatigue.Groups(iGroup).name)), globalOffset + [1 1]);
         case 2
            write_excel_cells( sheetObj, sprintf('%s Time Until Failure (s) for various S/N Curves', char(Fatigue.Groups(iGroup).name)), globalOffset + [1 1]);
         case 3
            write_excel_cells( sheetObj, sprintf('%s Lifetime Damage (-) without Goodman Correction for various S/N Curves', char(Fatigue.Groups(iGroup).name)), globalOffset + [1 1]);
         case 4
            write_excel_cells( sheetObj, sprintf('%s Time Until Failure (s) without Goodman Correction for various S/N Curves', char(Fatigue.Groups(iGroup).name)), globalOffset + [1 1]);
      end
      
      rangeObj = sheetObj.Range(convertR1C1toA1(globalOffset + [1 1], globalOffset + [1 2+nGroupChannels]) );
      rangeObj.Merge(1);
      rangeObj.Font.Bold = 1;
      rangeObj.HorizontalAlignment = -4108;
      rangeObj.VerticalAlignment = 2;
      rangeObj.RowHeight = 18;

      if (FileInfo.HaveUnits)
         unitsRow = 1;        
      else
         unitsRow = 0;
      end
      
      write_excel_cells( sheetObj, 'L_Ult', globalOffset + [3+unitsRow 1]);
      %rangeObj = sheetObj.Range(convertR1C1toA1([6 1], [6 1]) );
      % The line below should work but doesn't
      % rangeObj.Characters(2,3).Font.Superscript=1;
      
      rangeObj = sheetObj.Range(convertR1C1toA1(globalOffset + [3+unitsRow 1], globalOffset + [3+unitsRow 1]) );
      rangeObj.Font.Bold = 1;
      allSlopes = Fatigue.Groups(iGroup).allSNSlopes;
      nAllSlopes = length(allSlopes);
      write_excel_cells( sheetObj, 'm', globalOffset+[4+unitsRow 1]);
      rangeObj = sheetObj.Range(convertR1C1toA1(globalOffset + [4+unitsRow 1], globalOffset + [3+unitsRow+nAllSlopes 1]) );
      rangeObj.Merge(0);
      rangeObj.Font.Bold = 1;
      %rangeObj.Orientation = 90;
      rangeObj.ColumnWidth = 15;
      rangeObj.WrapText = 1;

      for i=1:nAllSlopes
         write_excel_cells( sheetObj, allSlopes(i), globalOffset + [3+unitsRow+i 2]);

      end
      rangeObj = sheetObj.Range(convertR1C1toA1(globalOffset + [4+unitsRow 2], globalOffset +[3+unitsRow+nAllSlopes 2]) );
      rangeObj.Font.Bold = 1;

      % Write out data input cells
      rangeObj = sheetObj.Range(convertR1C1toA1(globalOffset + [2 3], globalOffset + [2 2+nGroupChannels]) );
      rangeObj.Font.Bold = 1;
      rangeObj.ColumnWidth = 15;

      strData = cell(nFiles*nAllSlopes,nGroupChannels);

      for i=1:nGroupChannels
         write_excel_cells(sheetObj, FileInfo.Names{Fatigue.ChanInfo(Fatigue.Groups(iGroup).channelIndices(i)).Chan}, globalOffset + [2 2+i]);
         if (FileInfo.HaveUnits)
            write_excel_cells(sheetObj, FileInfo.Units{Fatigue.ChanInfo(Fatigue.Groups(iGroup).channelIndices(i)).Chan}, globalOffset + [3 2+i]);
            rangeObj = sheetObj.Range(convertR1C1toA1(globalOffset + [3 3], globalOffset + [3 2+nGroupChannels]) );
            rangeObj.Font.Bold = 1;
         end
         
         % write Lult
         rangeObj = write_excel_cells(sheetObj, sprintf( realFmt,Fatigue.ChanInfo(Fatigue.Groups(iGroup).channelIndices(i)).LUlt), globalOffset + [3+unitsRow 2+i]);
         format_number_cells(rangeObj,realFmt);
         
         % write DELs
         for j=1:Fatigue.ChanInfo(Fatigue.Groups(iGroup).channelIndices(i)).NSlopes
            rowOffset = find((Fatigue.Groups(iGroup).allSNSlopes == Fatigue.ChanInfo(Fatigue.Groups(iGroup).channelIndices(i)).SNslopes(j)) > 0);
            switch lifetimeType
                  case 1
                     value = Fatigue.Channel(Fatigue.Groups(iGroup).channelIndices(i)).lifetimeDamage(j);
                  case 2
                     value = Fatigue.Channel(Fatigue.Groups(iGroup).channelIndices(i)).timeUntilFailure(j);
                  case 3
                     value = Fatigue.Channel(Fatigue.Groups(iGroup).channelIndices(i)).lifetimeDamage_NoGoodman(j);
                  case 4
                     value = Fatigue.Channel(Fatigue.Groups(iGroup).channelIndices(i)).timeUntilFailure_NoGoodman(j);
            end
            
            %write_excel_cells( sheetObj, sprintf( realFmt, value ), localOffset );
            strData{rowOffset,i} = sprintf( realFmt, value );
         end
      end

      % Write the data as a single array
      rangeObj = write_excel_cells( sheetObj, strData, globalOffset + [4+unitsRow 3]);
      format_number_cells(rangeObj,realFmt);
      rangeObj = sheetObj.Range(convertR1C1toA1(globalOffset + [2 1], globalOffset + [3+unitsRow+nAllSlopes 2+nGroupChannels]));
      rangeObj.HorizontalAlignment = -4108;
      rangeObj.VerticalAlignment = 2;
      rangeObj.RowHeight = 18;
      rangeObj.Border.Value = 1;
      rangeObj.BorderAround(1,3);

      rangeObj = sheetObj.Range(convertR1C1toA1(globalOffset+[2+unitsRow, 1], globalOffset+[2+unitsRow, 2+nGroupChannels]) );
      rangeObj.Borders.Item(4).LineStyle = 9; % double line on bottom border
   
      globalOffset = globalOffset + [6+unitsRow+nAllSlopes 0];
   end