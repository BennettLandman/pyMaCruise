import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.List;
import edu.jhmi.rad.medic.structures.BinaryHeap4D;

/** Exercises the selected propagateThinning outer-loop condition. */
public final class ToadsThinningQueueRunReferenceHarness {
    private static String item(BinaryHeap4D queue) {
        return "["+Float.floatToIntBits(queue.getFirst())+","+queue.getFirstX()+","+
            queue.getFirstY()+","+queue.getFirstZ()+","+queue.getFirstK()+"]";
    }
    private static String run(String name,float limit,boolean suppressScheduling) throws Exception {
        Class<?> stepHarness=Class.forName("ToadsThinningQueueStepReferenceHarness");
        Class<?> stateType=Class.forName("ToadsThinningQueueStepReferenceHarness$State");
        Constructor<?> constructor=stateType.getDeclaredConstructor(String.class);
        constructor.setAccessible(true);
        Object state=constructor.newInstance("direct_allowed");
        Field queueField=stateType.getDeclaredField("queue");queueField.setAccessible(true);
        BinaryHeap4D queue=(BinaryHeap4D)queueField.get(state);
        if(suppressScheduling) {
            Field availableField=stateType.getDeclaredField("available");availableField.setAccessible(true);
            boolean[][][][] available=(boolean[][][][])availableField.get(state);
            for(int x=0;x<available.length;x++)for(int y=0;y<available[x].length;y++)
                for(int z=0;z<available[x][y].length;z++)for(int k=0;k<available[x][y][z].length;k++)
                    available[x][y][z][k]=false;
        }
        Method process=stepHarness.getDeclaredMethod("process",stateType);process.setAccessible(true);
        List<String> popped=new ArrayList<>(),snapshots=new ArrayList<>();float val=0.0f;
        while(queue.isNotEmpty() && val > -limit) {
            val=queue.getFirst();popped.add(item(queue));process.invoke(null,state);
            Field stepSeg=stateType.getDeclaredField("seg");stepSeg.setAccessible(true);
            Field stepSkel=stateType.getDeclaredField("skel");stepSkel.setAccessible(true);
            Method stepLabels=stepHarness.getDeclaredMethod("labels",int[][][].class);stepLabels.setAccessible(true);
            Method stepSkeleton=stepHarness.getDeclaredMethod("skeleton",int[][][].class);stepSkeleton.setAccessible(true);
            snapshots.add("{\"segmentation\":\""+stepLabels.invoke(null,stepSeg.get(state))+"\",\"skeleton\":\""+stepSkeleton.invoke(null,stepSkel.get(state))+"\"}");
        }
        Field segField=stateType.getDeclaredField("seg");segField.setAccessible(true);
        Field skelField=stateType.getDeclaredField("skel");skelField.setAccessible(true);
        Field avField=stateType.getDeclaredField("available");avField.setAccessible(true);
        Field volumeField=stateType.getDeclaredField("volume");volumeField.setAccessible(true);
        Method labels=stepHarness.getDeclaredMethod("labels",int[][][].class);labels.setAccessible(true);
        Method skeleton=stepHarness.getDeclaredMethod("skeleton",int[][][].class);skeleton.setAccessible(true);
        Method availability=stepHarness.getDeclaredMethod("availability",boolean[][][][].class);availability.setAccessible(true);
        Method floats=stepHarness.getDeclaredMethod("floats",float[].class);floats.setAccessible(true);
        Method queueJson=stepHarness.getDeclaredMethod("queue",BinaryHeap4D.class);queueJson.setAccessible(true);
        return "{\"name\":\""+name+"\",\"limit\":"+limit+",\"suppress_scheduling\":"+suppressScheduling+
            ",\"processed\":["+String.join(",",popped)+"],\"snapshots\":["+String.join(",",snapshots)+"],\"segmentation\":\""+
            labels.invoke(null,segField.get(state))+"\",\"skeleton\":\""+skeleton.invoke(null,skelField.get(state))+
            "\",\"available\":\""+availability.invoke(null,avField.get(state))+"\",\"volumes\":"+
            floats.invoke(null,volumeField.get(state))+",\"queue_after\":"+queueJson.invoke(null,queue)+"}";
    }
    public static void main(String[] args) throws Exception {
        ToadsHomeomorphicChangeReferenceHarness.initializeReferenceLut();
        System.out.print("{\"source\":\"VandyCRUISE MultipleToadSegmentation.propagateThinning outer queue loop\",\"source_sha256\":\"63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884\",\"cases\":[");
        System.out.print(run("drain_with_rescheduling",100.0f,false));
        System.out.print(","+run("cross_threshold_at_second_item",0.5f,true));
        System.out.print(","+run("top_priority_equals_threshold",0.25f,true));
        System.out.print(","+run("zero_limit_processes_nothing",0.0f,true));
        System.out.println("]}");
    }
}
