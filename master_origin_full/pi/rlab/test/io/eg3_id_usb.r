//
//
//

static(udevadm_propeties);
if (!exist(udevadm_propeties))
{
  udevadm_propeties = [ ...
      "BUSNUM", "DEVNAME", "DEVNUM", "DEVPATH", "DEVTYPE", "DRIVER", "ID_BUS", ...
      "ID_FOR_SEAT", "ID_MODEL", "ID_MODEL_ENC", "ID_MODEL_FROM_DATABASE", ...
      "ID_MODEL", "ID_PATH", "ID_PATH_TAG", "ID_REVISION", "ID_SERIAL", ...
      "ID_SERIAL_SHORT", "ID_USB_INTERFACES", "ID_VENDOR", "ID_VENDOR_ENC", ...
      "ID_VENDOR_FROM_DATABASE", "ID_VENDOR_ID", "MAJOR", "MINOR", ...
      "PRODUCT", "SUBSYSTEM", "TAGS", "TYPE", "USEC_INITIALIZED" ...
  ];
}

usbdevname2 = function(idlist)
{
  if (class(idlist)!="list")
  { return blank(0,0); }

  for (m in members(idlist))
  {
    i = find(strindex(toupper(m), udevadm_propeties));
    if(isempty(i))
    {
      printf("usbfname: Field '%m' is not recognized!\n", m);
      printf("usbfname: Available fields are (case is not important):\n");
      udevadm_propeties
      return blank(0,0);
    }
  }

  devs = reads("|find /sys/bus/usb/devices/usb*/ -name dev 2>/dev/null" );
  devs = rstrip(devs, "/dev");
  rval = blank(0,0);
  for (d in devs)
  {
    p = reads("|udevadm info -q property --export -p " + d + " 2>/dev/null");
    i = [];
    for (m in members(idlist))
    {
      i1 = find(strindex(p,toupper(m)));
      if (isempty(i1))
      { break; } // entry from idlist is not found, go to the next device
      for (i2 in i1)
      {
        if (strindex(p[i2], idlist.[m]))
        {
          i = [i2];
        else
          i = [];
        }
      }
      if (isempty(i))
      { break; }
    }
    if (isempty(i))
    { continue; }

    j = strindex(p,"DEVNAME=");
    if (isempty(j))
    { continue; }
    k = min(find(j));
    rval = [rval, p[k]];
  }

  if (exist(rval))
  {
    rval = lstrip(rval, "DEVNAME='");
    rval = rstrip(rval, "'");
  }

  // if the user has provided 'devname' or part of it
  // then use it to filter the output with
  m = members(idlist);
  i = find(strindex("DEVNAME", toupper(m)));
  if (!isempty(i))
  {
    j = find(strindex(rval, idlist.[m[i]]));
    if (isempty(j))
    {
      rval = blank(0,0);
    else
      rval = rval[j];
    }
  }

  return rval;
};


//
// find which port the device with properties below is attached to
//
attr=<<>>;
attr.devname = "tty";
attr.id_vendor_id = "0403";
attr.id_serial_short = "A600CJ5X";

tic();
r1 = usbdevname(attr)
toc()

tic();
r2 = usbdevname2(attr)
toc()
