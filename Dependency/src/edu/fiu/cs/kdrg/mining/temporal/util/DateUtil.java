package edu.fiu.cs.kdrg.mining.temporal.util;

import java.text.SimpleDateFormat;
import java.util.Date;

/**
 * 
 * @author Liang Tang
 * @date Dec 27, 2013 4:48:24 PM
 */
public class DateUtil {

  public final static SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd hh:mm:ss");

  public static String getDate() {
    return sdf.format(new Date());
  }

}
