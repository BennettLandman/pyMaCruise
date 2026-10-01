/** Pinned behavior of LesionToadPreprocess.updateTransformedTemplate. */
public final class LesionToadPreprocessUpdateTransformedTemplateReferenceHarness {
    private static final class State {
        byte[][][] topology;
        boolean transformed;
        void updateTransformedTemplate(byte[][][] tpl) {
            topology = tpl;
            transformed = true;
        }
    }

    public static void main(String[] args) {
        byte[][][] template = new byte[5][6][7];
        for (int x=0;x<5;x++) for(int y=0;y<6;y++) for(int z=0;z<7;z++)
            template[x][y][z] = (byte)(1 + ((x*11 + y*5 + z*3 + x*y) % 13));
        State state = new State();
        state.updateTransformedTemplate(template);
        boolean sameReference = state.topology == template;
        boolean sameContents = true;
        for(int x=0;x<5;x++) for(int y=0;y<6;y++) for(int z=0;z<7;z++)
            sameContents &= state.topology[x][y][z] == template[x][y][z];
        System.out.print("{\"source_method\":\"LesionToadPreprocess.updateTransformedTemplate\",\"shape\":[5,6,7],\"transformed\":"+state.transformed+",\"same_reference\":"+sameReference+",\"same_contents\":"+sameContents+",\"labels\":[");
        boolean first=true;
        for(int x=0;x<5;x++)for(int y=0;y<6;y++)for(int z=0;z<7;z++){
            if(!first)System.out.print(","); first=false;
            System.out.print(state.topology[x][y][z]);
        }
        System.out.println("]}");
    }
}
