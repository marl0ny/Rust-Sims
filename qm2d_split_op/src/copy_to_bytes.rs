/* The functions in this file transforms a variable or object 
into its bytes representation and dumps it to a byte array at the given
offset.
*/

pub fn copy_f32(bytes: &mut [u8], val: f32, offset: &mut usize) {
    let val_arr: [u8; 4] = val.to_ne_bytes();
    for i in 0..4 {
        bytes[*offset] = val_arr[i];
        *offset += 1;
    }
}

pub fn copy_i32(bytes: &mut [u8], val: i32, offset: &mut usize) {
    let val_arr: [u8; 4] = val.to_ne_bytes();
    for i in 0..4 {
        bytes[*offset] = val_arr[i];
        *offset += 1;
    }
}

pub fn copy_u32(bytes: &mut [u8], val: u32, offset: &mut usize) {
    let val_arr: [u8; 4] = val.to_ne_bytes();
    for i in 0..4 {
        bytes[*offset] = val_arr[i];
        *offset += 1;
    }
}

pub fn copy_u16(bytes: &mut [u8], val: u16, offset: &mut usize) {
    let val_arr: [u8; 2] = val.to_ne_bytes();
    for i in 0..2 {
        bytes[*offset] = val_arr[i];
        *offset += 1;
    }
}
