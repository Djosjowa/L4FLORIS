function write_file_statistics_as_excel(sheetObj, progName, iFile, FileInfo, Statistics, realFmt)

   % Get the date and time.

   DateTime = clock;

   globalOffset = [1 1];
   
   % Create table headers, etc.
   [~, matVer] = strtok(version,'(');
   write_excel_cells( sheetObj, sprintf( 'These statistics were generated by %s on %s at %02d:%02d:%02d by MATLAB %s.', progName, date, uint8( DateTime(4:6) ), matVer ), globalOffset );
   globalOffset = globalOffset + [1 0];
   
   write_excel_cells( sheetObj, sprintf( 'The analysis was based upon %d rows.', FileInfo.nSamples(iFile) ), globalOffset );
   globalOffset = globalOffset + [2 0];
   
   % Create table title using filename
   [~, name, ~] = fileparts(FileInfo.FileName{iFile});
   write_excel_cells( sheetObj, [name ' Statistics'], globalOffset );
   
   globalOffset = globalOffset + [1 0];
   
   % Create table Header
   write_excel_cells( sheetObj, 'Channel', globalOffset );
   headerOffset = globalOffset + [0 1];
   if (FileInfo.HaveUnits)
      write_excel_cells( sheetObj, 'Units', headerOffset );
      headerOffset = headerOffset + [0 1];
   end
   write_excel_cells( sheetObj, 'Minimum', headerOffset  );
   write_excel_cells( sheetObj, 'Mean', headerOffset + [0 1] );
   write_excel_cells( sheetObj, 'Maximum', headerOffset + [0 2] );
   write_excel_cells( sheetObj, 'StdDev', headerOffset + [0 3] );
   write_excel_cells( sheetObj, 'Skewness', headerOffset + [0 4] );
   write_excel_cells( sheetObj, 'Kurtosis', headerOffset + [0 5] );
   write_excel_cells( sheetObj, 'Range', headerOffset + [0 6] );
   
   
   % Create a master cell array for all statistics results
   strData = cell(FileInfo.nChannels,7);
   for iCh = 1:FileInfo.nChannels
      strData{iCh,1} = sprintf( realFmt, Statistics.Minima(iFile, iCh));
      strData{iCh,2} = sprintf( realFmt, Statistics.Means(iFile, iCh));
      strData{iCh,3} = sprintf( realFmt, Statistics.Maxima(iFile, iCh));
      strData{iCh,4} = sprintf( realFmt, Statistics.StdDevs(iFile, iCh));
      strData{iCh,5} = sprintf( realFmt, Statistics.Skews(iFile, iCh));
      strData{iCh,6} = sprintf( realFmt, Statistics.Kurtosis(iFile, iCh));
      strData{iCh,7} = sprintf( realFmt, Statistics.Range(iFile, iCh));
   end
   
   % Write channel data
   for iCh = 1:FileInfo.nChannels
      % Channel name
      write_excel_cells( sheetObj, FileInfo.Names{iCh}, globalOffset + [iCh 0] );
         
            
      if (FileInfo.HaveUnits)
         % Units
         write_excel_cells( sheetObj, FileInfo.Units{iCh}, globalOffset + [iCh 1] );
      end
%       % Min
%        write_excel_cells( sheetObj, sprintf( realFmt, Statistics.Minima(iFile, iCh)), headerOffset + [iCh 0] );
%       % Mean
%       write_excel_cells( sheetObj, sprintf( realFmt, Statistics.Means(iFile, iCh)), headerOffset + [iCh 1] );
%       % Max
%       write_excel_cells( sheetObj, sprintf( realFmt, Statistics.Maxima(iFile, iCh)), headerOffset + [iCh 2] );
%       % StdDev
%       write_excel_cells( sheetObj, sprintf( realFmt, Statistics.StdDevs(iFile, iCh)), headerOffset + [iCh 3] );
%       % Skew
%       write_excel_cells( sheetObj, sprintf( realFmt, Statistics.Skews(iFile, iCh)), headerOffset + [iCh 4] );
%        % Kurtosis
%       write_excel_cells( sheetObj, sprintf( realFmt, Statistics.Kurtosis(iFile,iCh)), headerOffset + [iCh 5] );
%       % Range
%       write_excel_cells( sheetObj, sprintf( realFmt, Statistics.Range(iFile, iCh)), headerOffset + [iCh 6] );
   end
   write_excel_cells( sheetObj, strData, headerOffset + [1 0]);
   
   % Format table cells
   rangeObj = sheetObj.Range( convertR1C1toA1(headerOffset, headerOffset + [FileInfo.nChannels 6]) );
   format_number_cells(rangeObj, realFmt);
   
   % Format Column headers and table title
   rangeObj = sheetObj.Range( convertR1C1toA1([4 1], [4 headerOffset(2)+6]) );
   rangeObj.Font.Bold = 1;
   rangeObj.Merge(0);
   rangeObj.HorizontalAlignment = -4108;
   rangeObj = sheetObj.Range( convertR1C1toA1([5 1], [5 headerOffset(2)+6]) );
   rangeObj.Font.Bold = 1;
   
   % Format table cells
   rangeObj = sheetObj.Range(convertR1C1toA1(globalOffset, globalOffset + [FileInfo.nChannels headerOffset(2)+5]) );
   rangeObj.VerticalAlignment = 2;
   rangeObj.RowHeight = 18;
   rangeObj.Columns.AutoFit();
   rangeObj.HorizontalAlignment = -4108;
   rangeObj.Border.Value = 1;
   rangeObj.BorderAround(1,3);   
end