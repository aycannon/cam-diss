# cam-diss
Quantile Analysis Approach to the Forward premium Bias


---
#### Notes
The data is pulled from the Refinitiv LSEG Workspace
- The data is stored in the `data` folder
Function form for the data pull in Excel is:

`=@RDP.Data("GBP1MV=;EUR1MV=...", "TR.BIDPRICE;TR.ASKPRICE;TR.MIDPRICE;TR.MIDPRICE.DATE", "SDate=2000-01-01 Frq=D EDate=0 CH=Fd RH=IN")`
EDate = 0 just means the end date is today. Getting the date from this TR.MIDPRICE.DATE function.

[Loss Function Animation](loss_animation.gif)