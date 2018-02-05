# Image-processing
環境 visual studio 2015 + OpenCV 2.4

操作功能
1.把filename改成要開啟的圖片名稱，stainedcell_small或者stainedcell_large
2.調整Threshold，small時用195，large時用190
3.結果會放在Project資料夾下的filename資料夾，裡面有txt的結果檔與輸出圖

GUI功能介紹
1.對著原圖點擊滑鼠左鍵，會顯示出最近的點擊處的細胞，其相關資訊會列在左上角
2.原圖上方會有Trackbar，拉動可以調整放大倍率，會把點擊目標細胞作放大，倍率範圍為1~10
3.若選到的細胞是想要紀錄的，滑鼠點選右鍵就可以儲存
